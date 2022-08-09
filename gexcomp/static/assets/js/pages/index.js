

var mainChart = null;
var mainChartWinStart = 0
var mainChartWinEnd = 100

function createMainChart(data) {

    mainChart = Highcharts.stockChart('main-chart', {

        rangeSelector: {
            selected: 1
        },

        xAxis: {
            events: {
                setExtremes: (event) => {
                    mainChartWinStart = Math.round(event.min);
                    mainChartWinEnd = Math.round(event.max);

                    $('#win-len').text(mainChartWinEnd - mainChartWinStart + 1);
                    //console.log(mainChartWinStart + " - " + mainChartWinEnd);
                    //event.preventDefault();
                }
            },
        },

        yAxis: {
            labels: {
                formatter: function () {
                    return this.value;
                }
            },
            plotLines: [{
                value: 0,
                width: 2,
                color: 'silver'
            }]
        },

        plotOptions: {
            series: {
                //compare: 'percent',
                showInNavigator: true
            }
        },

        tooltip: {
            pointFormat: '<span style="color:{series.color}">{series.name}</span>: <b>{point.y}</b><br/>',
            valueDecimals: 2,
            split: true
        },

        series: data['chart-series'],

        navigator: {
            series: {
                data: data['nav-series']
            }
        },
    });
}

function getPointCategoryName(point, dimension) {
    var series = point.series,
        isY = dimension === 'y',
        axis = series[isY ? 'yAxis' : 'xAxis'];
    return axis.categories[point[isY ? 'y' : 'x']];
}

function createWinHeatmap(data) {

    Highcharts.chart('gex-err-heatmap', {

        chart: {
            type: 'heatmap',
            marginTop: 40,
            marginBottom: 80,
            plotBorderWidth: 1,
            height: 200
        },


        title: {
            text: ''
        },

        xAxis: {
            categories: data['x-axis']
        },

        yAxis: {
            categories: ['Prediction Error'],
            title: 'Prediction Error',
            reversed: true
        },

        accessibility: {
            point: {
                descriptionFormatter: function (point) {
                    var ix = point.index + 1,
                        geneName = getPointCategoryName(point, 'x'),
                        yName = getPointCategoryName(point, 'y'),
                        val = point.value;
                    return geneName + ' pred. err = ~' + val + '%';
                }
            }
        },

        colorAxis: {
            min: 0,
            minColor: '#FFFFFF',
            maxColor: Highcharts.getOptions().colors[0]
        },

        legend: {
            align: 'right',
            layout: 'vertical',
            margin: 0,
            verticalAlign: 'top',
            y: 25,
            symbolHeight: 280
        },

        tooltip: {
            formatter: function () {
                return '<b>' + getPointCategoryName(this.point, 'x') + '</b><br><b>Pred. Err: </b> ~' +
                    this.point.value + '%';
            }
        },

        series: [{
            name: 'pred-err',
            borderWidth: 1,
            data: data['pred-err-series'], //[[0, 0, 10], [0, 1, 19], [0, 2, 8], [0, 3, 24], [0, 4, 67], [1, 0, 92], [1, 1, 58], [1, 2, 78], [1, 3, 117], [1, 4, 48], [2, 0, 35], [2, 1, 15], [2, 2, 123], [2, 3, 64], [2, 4, 52], [3, 0, 72], [3, 1, 132], [3, 2, 114], [3, 3, 19], [3, 4, 16], [4, 0, 38], [4, 1, 5], [4, 2, 8], [4, 3, 117], [4, 4, 115], [5, 0, 88], [5, 1, 32], [5, 2, 12], [5, 3, 6], [5, 4, 120], [6, 0, 13], [6, 1, 44], [6, 2, 88], [6, 3, 98], [6, 4, 96], [7, 0, 31], [7, 1, 1], [7, 2, 82], [7, 3, 32], [7, 4, 30], [8, 0, 85], [8, 1, 97], [8, 2, 123], [8, 3, 64], [8, 4, 84], [9, 0, 47], [9, 1, 114], [9, 2, 31], [9, 3, 48], [9, 4, 91]],
            dataLabels: {
                enabled: true,
                color: '#000000'
            }
        }],

        responsive: {
            rules: [{
                condition: {
                    maxWidth: 500
                },
                chartOptions: {
                    yAxis: {
                        labels: {
                            formatter: function () {
                                return this.value.charAt(0);
                            }
                        }
                    }
                }
            }]
        }

    });


}


//This function is disabled
function createWordCloud(containerID, data) {

    return;

    const text = 'Lorem ipsum dolor sit amet, consectetur adipiscing elit. Aenean bibendum erat ac justo sollicitudin, quis lacinia ligula fringilla. Pellentesque hendrerit, nisi vitae posuere condimentum, lectus urna accumsan libero, rutrum commodo mi lacus pretium erat. Phasellus pretium ultrices mi sed semper. Praesent ut tristique magna. Donec nisl tellus, sagittis ut tempus sit amet, consectetur eget erat. Sed ornare gravida lacinia. Curabitur iaculis metus purus, eget pretium est laoreet ut. Quisque tristique augue ac eros malesuada, vitae facilisis mauris sollicitudin. Mauris ac molestie nulla, vitae facilisis quam. Curabitur placerat ornare sem, in mattis purus posuere eget. Praesent non condimentum odio. Nunc aliquet, odio nec auctor congue, sapien justo dictum massa, nec fermentum massa sapien non tellus. Praesent luctus eros et nunc pretium hendrerit. In consequat et eros nec interdum. Ut neque dui, maximus id elit ac, consequat pretium tellus. Nullam vel accumsan lorem.',
    lines = text.split(/[,\. ]+/g),
    dataTmp = lines.reduce((arr, word) => {
    let obj = Highcharts.find(arr, obj => obj.name === word);
    if (obj) {
      obj.weight += 1;
    } else {
      obj = {
        name: word,
        weight: 1
      };
      arr.push(obj);
    }
    return arr;
    }, []);

    Highcharts.chart(containerID, {
        accessibility: {
            screenReaderSection: {
                beforeChartFormat: '<h5>{chartTitle}</h5>' +
                    '<div>{chartSubtitle}</div>' +
                    '<div>{chartLongdesc}</div>' +
                    '<div>{viewTableButton}</div>'
            }
        },
        series: [{
            type: 'wordcloud',
            data: dataTmp,
            name: 'go-terms'
        }],
        title: {
            text: ''
        }
    });

    console.log("hi");
}

function createEnrichmentBarChart(containerID, data, regulation) {//regulation can be up or down

    if (["up-regulated", "down-regulated"].includes(regulation) == false) {
        alert("Error: wrong regulation type for encrichment bar chart. Must be 'up' or 'down'.")
        return;
    }

    let TEST_MODE = true;
    //categories
    var categories = (TEST_MODE) ? ["(term1)", "(t2)", "(t3)", "(t4)"] : data[regulation]["terms"];

    Highcharts.chart(containerID, {
        chart: {
            type: 'bar'
        },
        title: {
            text: ''
        },
        subtitle: {
            text: (TEST_MODE) ? 'Barchart some DB' : data["database"]
        },
        xAxis: [{
            categories: categories,
            reversed: false,
            labels: {
                step: 1
            }
        }],
        yAxis: {
            title: {
                text: null
            },
            labels: {
                formatter: function () {
                    return '';//parseFloat(this.value).toFixed(3);
                }
            }
        },

        plotOptions: {
            series: {
                stacking: 'normal'
            }
        },

        tooltip: {
        formatter: function () {
            return '<b>' + this.series.name + ', Term: ' + this.point.category + '</b><br/>' +
            '-log10 P-val: ' +  -1 * Math.log10(parseFloat(this.point.y));
            }
        },

        series: [
            {
                name: 'Terms from ' + regulation + ' genes',
                data: (TEST_MODE) ? [1,1.5,1.2,1.4] : data[regulation]["pvals"],
                cursor: 'pointer',
                color: (regulation == "up-regulated") ? '#F8766D' : '#619CFF',
                point: {
                    events: {
                        click: function () {
                            //TODO: add functionality to run whole profile enrichment analysis for selected term
                            alert('Category: ' + this.category.substring(this.category.indexOf('(') + 1, this.category.indexOf(')')));
                        }
                    }
                }
            }
        ]
    });

}

function createColChart(containerID, data, yLabelsEnabled = false) {//regulation can be up or down

    Highcharts.chart(containerID, {
        chart: {
            type: 'column'
        },
        title: {
            text: 'Count of Chr-DEGs in the Selected Window'
        },
        subtitle: {
            text: ''
        },
        xAxis: {
            categories: data["chromosomes"]
        },
        yAxis: {
            title: {
                text: ''
            },
            labels: {
                enabled: yLabelsEnabled
            }

        },
        plotOptions: {
            series: {
                stacking: 'normal'
            }
        },

        tooltip: {
        formatter: function () {
                return '<b>Chr: ' + this.point.category + '</b><br/>' + this.point.y + ' DEG(s)';
            }
        },

        series: [
            {
                name: 'Chr-DEGs Count ',
                colorByPoint: true,
                data: data['series']
            }
        ]
    });

}

function displayWinOnNavigator(winStart, winEnd) {

    //scroll to page top
    $('html,body').animate({ scrollTop: 0 }, 'slow', function () {
        //move navigaor to selected win
        mainChart.xAxis[0].setExtremes(winStart, winEnd);
    });



}


function startWinSelectedBioAnalysis() {
    let runID = "1";

    //win pred. err. heatmap
    //creat window heatmap
    json_window_heatmap_url = `/json_win_heatmap/${runID}/${mainChartWinStart}/${mainChartWinEnd}`
    $.get(json_window_heatmap_url, function(response) {
    })
        .done(function(response) {
            if (response) {

                let data = response;
                console.log(data);

                //change values of pred-err-mean and pred-err-sum spans
                let predErrMean = data["pred-err-mean"];
                let predErrSum = data["pred-err-sum"];
                $('#pred-err-mean').text(predErrMean);
                $('#pred-err-sum').text(predErrSum);

                createWinHeatmap(data);

                //show win heatmap div
                $('#win-heatmap').removeClass('d-none');

                //scroll to heatmap
                $('html, body').animate({ scrollTop:$('#win-heatmap').position().top }, 'slow', function () {});
            }
        })
        .fail(function() {
            alert( "Error in fetching heatmap data. Please try again later." );
        })

    console.log("bio win selected analysis starting..");


    json_win_selected_bio_analysis = `/json_win_selected_bio_analysis/${runID}/${mainChartWinStart}/${mainChartWinEnd}`;
    $.get(json_win_selected_bio_analysis, function(response) {
    })
      .done(function(response) {
          if (response) {
              let data = response;
              console.log(data);

              //create selected win GO word cloud & barchart
              let dataGOBarChart = data['data-go']['bar-chart']; //will contain data for up and down regulated
              createEnrichmentBarChart('win-go-up-barchart', dataGOBarChart, 'up-regulated');
              createEnrichmentBarChart('win-go-down-barchart', dataGOBarChart, 'down-regulated');
              //createWordCloud('win-go-word-cloud', dataGOWordCloud);

              //create selected win KEGG analysis word cloud & barchart
              let dataKEGGChartBar = data['data-kegg']['bar-chart']; //will contain data for up and down regulated
              //let dataKEGGWordCloud = '';
              createEnrichmentBarChart('win-kegg-up-barchart', dataKEGGChartBar, 'up-regulated');
              createEnrichmentBarChart('win-kegg-down-barchart', dataKEGGChartBar, 'down-regulated');
              //createWordCloud('win-kegg-word-cloud', dataKEGGWordCloud);

              let dataTblIncludedGenes = data['tbl-included-genes'];
              $("#tbl-win-included-genes > tbody").html(dataTblIncludedGenes);

              let dataTblTopDEGs = data['tbl-top-degs'];
              $('#tbl-top-win-degs > tbody').html(dataTblTopDEGs);

              let dataDEAChrGeneCnt = data['data-dea-chr-gene-cnt'];
              createColChart('win-chr-degs-counts-colchart', dataDEAChrGeneCnt);

              //show win-info sidebar
              $('#win-info').removeClass('d-none');

              console.log("bio win selected analysis end.");
          }
      })
      .fail(function() {
            alert( "Error in starting the biological analysis. Please try again later." );
      })

}

function startWinLenBioAnalysis() {
    let runID = "1";
    let n_wins = parseInt($('#txt-n-wins').val());
    let win_len = mainChartWinEnd - mainChartWinStart + 1;
    let win_step = 5;

    console.log("bio win len analysis starting..");

    json_win_len_bio_analysis = `/json_win_len_bio_analysis/${runID}/${n_wins}/${win_len}/${win_step}`;
    $.get(json_win_len_bio_analysis, function(response) {

    })
      .done(function(response) {
          if (response) {
              let data = response;
              console.log(data);

              topWinsTblHTML = data["top-wins-tbl-html"];
              $("#tbl-win-len > tbody").html(topWinsTblHTML);

              //show top-wins-info div
              $('#top-wins-info').removeClass("d-none");

              console.log("bio win len analysis end.");
          }

      })
      .fail(function() {
            alert( "Error in starting the biological analysis. Please try again later." );
      })
}


$(document).ready(function() {
    //creat main-chart
    json_main_chart_url = "/json_main_chart"
    $.get(json_main_chart_url, function(response) {
    })
      .done(function(response) {
          if (response) {
              let data = response;
              console.log(data);
              createMainChart(data);
          }
      })
      .fail(function() {
            alert( "Error in fetching main chart data. Please try again later." );
      })



    //bind events

    $('#btn-win-selected-bio-analysis').click(function() {
        startWinSelectedBioAnalysis();
    });

    $('#btn-win-len-bio-analysis').click(function() {
        startWinLenBioAnalysis();
    });


});

