# -*- encoding: utf-8 -*-
"""
Copyright (c) 2019 - present AppSeed.us
"""

from django import template
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpResponseRedirect, JsonResponse
from django.template import loader
from django.urls import reverse
import config.settings.base as settings

import pandas as pd
import numpy as np
import random

from os import path
#to run bash scripts
import subprocess
from subprocess import Popen, PIPE

#EnrichR wrapper
import gseapy as gp

from enum import Enum
class SortGenesResult(Enum):
    lstGeneNames = 1
    dfGenesDetails = 2


hgnc_complete_set_path = path.join(settings.DATA_DIR, 'hgnc_complete_set_2022-07-01.csv')

def util_is_number(value):
    if type(value) == int or type(value) == float:
        return value
    else:
        raise Exception("The value you passed ({0}) is not a number.".format(value))

#TODO: test this function
#sorts the  provided Ensemble Gene IDs according to HGNC complete set. Genes that are not found are moved to the end of the set.
def util_sort_genes(gene_ens_ids, sort_genes_result):
    # HGNC complete gene set sorted by chr -> arm (p -> q) -> pos
    df_complete = pd.read_csv(hgnc_complete_set_path, index_col=False)
    df_complete = df_complete[["Sybmol", "chr", "arm", "pos", "Ensembl gene ID"]]
    # zip 2 cols as dict. ex ENSG ID -> 3 (chrom. 3)
    complete_chr_gene = dict(zip(df_complete["Ensembl gene ID"], df_complete["chr"]))

    sorted_genes = []
    unfound_genes = []
    for gene_id in gene_ens_ids:
        try:
            result = df_complete[df_complete["Ensemble gene ID"]].loc[0]
            #remove old order list
            gene_ens_ids.remove(gene_id)
            sorted_genes.append([gene_id, result["Approved symbol"], result["chr"], result["arm"], result["pos"]])
        except: #ID not found in complete set
            unfound_genes.append([gene_id, "NA", "NA", "NA", "NA"])

    #put unfound at the end of sorted genes
    sorted_genes.extend(unfound_genes)

    df_sorted = pd.DataFrame(unfound_genes, columns=["Ensemble gene ID", "symbol", "chr", "arm", "pos"])

    if sort_genes_result == sort_genes_result.dfGenesDetails:
        return df_sorted
    elif sort_genes_result == sort_genes_result.lstGeneNames:
        return df_sorted["Ensemble gene ID"].values.tolist()


#TODO: see if you want to integrate this function with sorting function , since most of analysis use the symbol not the ENSG ID
#gets genes symbols from ENSG IDs
#a non-found gene will be removed from passed list
def util_query_genes_symbols(genes_ensg_ids):
    # HGNC complete gene set sorted by chr -> arm (p -> q) -> pos
    df_complete = pd.read_csv(hgnc_complete_set_path, index_col=False)
    df_complete = df_complete[["Approved symbol", "Ensembl gene ID"]]

    result = []

    for gene_id in genes_ensg_ids:
        try:
            symbol = df_complete[df_complete["Ensembl gene ID"] == gene_id].iloc[0,0]
            result.append([gene_id, symbol]) #we want to return only genes with symbols
        except:
            pass

    df_result = pd.DataFrame(result, columns=["Ensemble gene ID", "Approved symbol"])

    return df_result

class EnrichRAnalysisInputType(Enum):
    EnsemblGeneIDs = 1
    geneSymbols = 2

#runs EnrichR GO & KEGG analyis on list of gene symbols and returns 2 DFs, one for each of analysis with cols :(Gene_set, Term, Overlap, P-value, Adjusted P-value, Old P-value, Old Adjusted P-value, Odds Ratio, Combined Score, Genes)
def util_run_enrichr_analysis(gene_list, enrichr_analysis_input_type):

    #select db from default is human
    gp.get_library_name()

    if enrichr_analysis_input_type == EnrichRAnalysisInputType.EnsemblGeneIDs: #convert to gene symbols
        df_query_res = util_query_genes_symbols(gene_list)
        gene_list = df_query_res["Approved symbol"].values.tolist()
        print("############OOOOOO#############", gene_list)

    enr = gp.enrichr(gene_list= gene_list, #expects string input or file
                     gene_sets=['GO_Biological_Process_2018', 'KEGG_2016'],
                     organism='human',  # don't forget to set organism to the one you desired! e.g. Yeast
                     outdir=None,  # don't write to disk
                     )

    df_res = enr.results
    df_res_go = df_res[df_res["Gene_set"] == "GO_Biological_Process_2018"]
    df_res_kegg = df_res[df_res["Gene_set"] == "KEGG_2016"]

    return df_res_go, df_res_kegg

#@login_required(login_url="/login/")
def index(request):
    context = {'segment': 'index'}

    html_template = loader.get_template('home/index.html')
    return HttpResponse(html_template.render(context, request))

def json_main_chart(request):

    run_id = "1"

    n_genes = 4000#3000

    # we will get the aggregated data as group mean because we want to show only one curve for each plot
    df = pd.read_csv(path.join(settings.RUNS_DIR, run_id, "data_aggregated.csv"), index_col=0)

    gex_1 = df.iloc[0,:n_genes].values
    gex_2 = df.iloc[1,:n_genes].values

    imputation_err = np.abs(gex_1 - gex_2)

    #HGNC complete gene set sorted by chr -> arm (p -> q) -> pos
    # chr-gene
    df_complete = pd.read_csv(hgnc_complete_set_path, index_col=False)
    df_complete = df_complete[["chr", "Ensembl gene ID"]]
    #zip 2 cols as dict. ex ENSG ID -> 3 (chrom. 3)
    complete_chr_gene = dict(zip(df_complete["Ensembl gene ID"], df_complete["chr"]))

    upload_genes = list(df.columns)
    upload_chr_gene = {}
    chrs = []
    #init and fill with leading zeros. e.g: (01) instead of (1). Thihs is to match .csv file content
    for i in range(1, 23, 1):
        chrs.append(str(i).zfill(2))
    chrs = chrs + ["X", "Y"]

    #init all with 0
    for chr in chrs:
        upload_chr_gene[chr] = 0

    #count
    for gene in upload_genes:
        if gene in complete_chr_gene:
            chr = complete_chr_gene[gene]
            upload_chr_gene[chr] += 1

    context = {
        'chart-series': [
            {'name': 'GEX1', 'data': gex_1.tolist()},
            {'name': 'GeX2', 'data': gex_2.tolist()}
        ],
        'nav-series': imputation_err.tolist(),
        'chr-gene-count': upload_chr_gene
    }

    return JsonResponse(context)

def json_win_heatmap(request, run_id, win_start, win_end):

    win_start = int(util_is_number(win_start))
    win_end = int(int(util_is_number(win_end)))
    n_genes = win_end - win_start

    df = pd.read_csv(path.join(settings.RUNS_DIR, run_id, "data_aggregated.csv"), index_col=0)

    # select window
    df_win = df.iloc[:, win_start:win_end + 1]

    gex_1 = df_win.iloc[0,:].values
    gex_2 = df_win.iloc[1,:].values

    pred_err = np.abs(gex_1 - gex_2)

    pred_err_mean = round(np.mean(pred_err), 3)
    pred_err_sum = round(np.sum(pred_err), 3)

    #normalize [0,1]
    pred_err_norm = list((pred_err - np.min(pred_err)) / (np.max(pred_err) - np.min(pred_err)))

    x_axis = list(df.columns[win_start:win_end])

    #show as [0,100]
    pred_err_series = [[i, 0, int(pred_err_norm[i] * 100)] for i in range(len(pred_err_norm))]

    context = {
        'x-axis': x_axis,
        'pred-err-series': pred_err_series,
        'pred-err-mean': pred_err_mean,
        'pred-err-sum': pred_err_sum
    }

    return JsonResponse(context)

def json_win_selected_bio_analysis(request, run_id, win_start, win_end):

    win_start = int(util_is_number(win_start))
    win_end = int(int(util_is_number(win_end)))
    n_genes = win_end - win_start

    df = pd.read_csv(path.join(settings.RUNS_DIR, run_id, "data.csv"), index_col=0)
    genes = list(df.columns)[win_start:win_end + 1]


    #select window
    df_win = df.iloc[:, win_start:win_end + 1]

    #transpose
    df_win = df_win.transpose()

    #change index name from sample to gene after transpose
    df_win.index.names= ["gene"]

    #round all values to int. otherwise DESeq2 will not work
    df_win[df_win.columns] = df_win[df_win.columns].astype(int)

    #save norm counts
    df_win_path = path.join(settings.RUNS_DIR, run_id, "data_norm_counts_win_{0}_to_{1}.csv".format(win_start, win_end))
    df_win.to_csv(df_win_path, index=True)

    #subprocess.run(
    #    ["conda", "run", "-n", "r-env", "Rscript", path.join(settings.SCRIPTS_DIR, "win_selected_bio_analysis.R"),
    #     run_id, str(win_start), str(win_end)])

    # p = subprocess.Popen(["conda" , "run", "-n", "r-env", "Rscript", path.join(settings.SCRIPTS_DIR, "win_selected_bio_analysis.R"), run_id, str(win_start), str(win_end)]
    #                      , stdout=PIPE, stderr=PIPE)
    # stdout, stderr = p.communicate()

    # if stderr:
    #     raise Exception("Error " + str(stderr))

    #GO & KEGG
    df_go, df_kegg = util_run_enrichr_analysis(genes, EnrichRAnalysisInputType.EnsemblGeneIDs)

    #get top GO terms list created by R script
    #.....
    # word_cloud = [
    #     {
    #         "name": "Signal Transduction",
    #         "weight": 5
    #     },
    #     {
    #         "name": "Phosphorylation",
    #         "weight": 4
    #     }
    # ]

    # HTML fo tbl includes ENSG ID and symbol
    tbl_included_genes = ""
    gene_symbols = list(util_query_genes_symbols(genes).to_numpy()) #note that un-found genes will be removed from the set. result shape: [["Ensembl gene ID", "Approved symbol"], [], ..]
    for entry in gene_symbols:
        tbl_included_genes += "<tr>"
        tbl_included_genes += "<td class=\"text-gray-900\">" + entry[0] + "</td>"
        tbl_included_genes += "<td class=\"text-gray-500\">" + entry[1] + "</td>"
        tbl_included_genes += "<td class=\"text-gray-500\">" + "0" + "</td>" #STD .. for later #TODO
        tbl_included_genes += "</tr>"

    go_dict = {
        'bar-chart': {
            'database': df_go.iloc[0]["Gene_set"],
            'terms': df_go.iloc[0:10]["Term"].values.tolist(),
            'pvals': df_go.iloc[0:10]["P-value"].values.tolist()
        }
    }

    kegg_dict = {
        'bar-chart': {
            'database': df_kegg.iloc[0]["Gene_set"],
            'terms': df_kegg.iloc[0:10]["Term"].values.tolist(),
            'pvals': df_kegg.iloc[0:10]["P-value"].values.tolist()
        }
    }

    context = {
        'segment': 'win-selected-bio-analysis',
        'run-id': run_id,
        'win-start': win_start,
        'win-end': win_end,
        'tbl-included-genes': tbl_included_genes,
        'data-go': go_dict,
        'data-kegg': kegg_dict
        #'std-out': str(stdout)
    }

    return JsonResponse(context)


def json_win_len_bio_analysis(request, run_id, n_wins, win_len, win_step):

    n_wins = int(util_is_number(n_wins))
    win_len = int(util_is_number(win_len))
    win_step = int(util_is_number(win_step))

    df = pd.read_csv(path.join(settings.RUNS_DIR, run_id, "data_aggregated.csv"), index_col=0)

    n_genes = len(df.columns)

    n_wins = n_wins if n_genes % win_step == 0 else n_wins + 1
    print(n_genes, n_wins)

    #win-pred-err dict
    #stores key: win_start -> value: sum of pred. err.
    win_pred_err = {}

    for i in range(int(n_genes / win_step) if n_genes % win_step == 0 else int(n_genes / win_step) + 1   ):
        win_start = i * win_step
        win_end = win_start + win_step

        #handling last window
        if (win_end + 1) > n_genes:
            win_end = n_genes - 1
            win_start = win_end - win_step

        exp_orig = df.iloc[0, win_start:win_end].values
        exp_pred = df.iloc[1, win_start:win_end].values

        exp_err = np.sum(np.abs(exp_orig - exp_pred))
        win_pred_err[win_start] = int(exp_err)

    #let's do this only once
    win_pred_err_keys = list(win_pred_err.keys())
    win_pred_err_values = list(win_pred_err.values())

    # let's find the n_wins with highest pred.err
    # using sorted() + lambda + list slicing
    indices = sorted(range(len(win_pred_err_values)), key=lambda sub: win_pred_err_values[sub])[-n_wins:]
    indices.reverse() # .reverse() so we start with highest pred. err. first
    top_wins = []
    for idx in indices:
        top_wins.append(
            {
                "win_start": win_pred_err_keys[idx],
                "win_end": win_pred_err_keys[idx] + win_len,
                "pred_err": win_pred_err_values[idx]
            }
        )

    # generate html table
    tbl_html = ""
    for i in range(len(top_wins)):
        win = top_wins[i]
        tbl_html += "<tr>"

        #Nr.
        tbl_html += "<td class =\"text-gray-900\" scope=\"row\">"
        tbl_html += str(i + 1)
        tbl_html += "</td>"

        # SUM OF PRED. ERROR
        tbl_html += "<td class =\"text-gray-500\" scope=\"row\">"
        tbl_html += str(win["pred_err"])
        tbl_html += "</td>"

        # START POS.
        tbl_html += "<td class =\"text-gray-500\" scope=\"row\">"
        tbl_html += str(win["win_start"])
        tbl_html += "</td>"

        # END POS.
        tbl_html += "<td class =\"text-gray-500\" scope=\"row\">"
        tbl_html += str(win["win_end"])
        tbl_html += "</td>"

        # DETAILS
        tbl_html += "<td class=\"text-gray-500\" scope=\"row\">"
        tbl_html += "<a class=\"text-secondary lnk-win-info\" href='javascript:;' data-win-start='{0}' data-win-end='{1}' onclick=\"displayWinInfo({0}, {1})\">(click)</a>".format(win["win_start"], win["win_end"])
        tbl_html += "</td>"

        tbl_html += "</tr>"


    # p = subprocess.Popen(["conda" , "run" "-n r-env", "Rscript", path.join(settings.SCRIPTS_DIR, "win_len_bio_analysis.R"), run_id, str(win_len)]
    #                      , stdout=PIPE, stderr=PIPE)
    # stdout, stderr = p.communicate()
    #
    # if stderr:
    #     raise Exception("Error " + str(stderr))

    context = {
        'segment': 'win-len',
        'run-id': run_id,
        'win-len': win_len,
        'top-wins-tbl-html': tbl_html
    }

    return JsonResponse(context)

#@login_required(login_url="/login/")
def pages(request):
    context = {}
    # All resource paths end in .html.
    # Pick out the html file name from the url. And load that template.
    try:

        load_template = request.path.split('/')[-1]

        if load_template == 'admin':
            return HttpResponseRedirect(reverse('admin:index'))
        context['segment'] = load_template

        html_template = loader.get_template('home/' + load_template)
        return HttpResponse(html_template.render(context, request))

    except template.TemplateDoesNotExist:

        html_template = loader.get_template('home/page-404.html')
        return HttpResponse(html_template.render(context, request))

    except:
        html_template = loader.get_template('home/page-500.html')
        return HttpResponse(html_template.render(context, request))
