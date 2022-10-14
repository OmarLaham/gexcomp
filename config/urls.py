from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.urls import include, path, re_path
from django.views import defaults as default_views
from django.views.generic import TemplateView


from gexcomp import views


urlpatterns = [
    path("", views.index, name='home'),
    path("run/<str:run_id>", views.run, name='run'),
    path("json_main_chart/<str:run_id>", views.json_main_chart, name='json_main_chart'),
    path("json_win_heatmap/<str:run_id>/<int:win_start>/<int:win_end>", views.json_win_heatmap, name='json_window_heatmap'),
    path("json_win_selected_bio_analysis/<str:run_id>/<int:win_start>/<int:win_end>", views.json_win_selected_bio_analysis, name='json_win_selected_bio_analysis'),
    path("json_win_len_bio_analysis/<str:run_id>/<int:n_wins>/<int:win_len>/<int:win_step>", views.json_win_len_bio_analysis, name='json_win_len_bio_analysis'),
    # Matches any html file
    re_path(r'^.*\.*', views.pages, name='pages'),
    #path(
    #    "about/", TemplateView.as_view(template_name="pages/about.html"), name="about"
    #),
    # Django Admin, use {% url 'admin:index' %}
    path(settings.ADMIN_URL, admin.site.urls),
    # User management
    path("users/", include("gexcomp.users.urls", namespace="users")),
    path("accounts/", include("allauth.urls")),
    # Your stuff: custom urls includes go here
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)


if settings.DEBUG:
    # This allows the error pages to be debugged during development, just visit
    # these url in browser to see how these error pages look like.
    urlpatterns += [
        path(
            "400/",
            default_views.bad_request,
            kwargs={"exception": Exception("Bad Request!")},
        ),
        path(
            "403/",
            default_views.permission_denied,
            kwargs={"exception": Exception("Permission Denied")},
        ),
        path(
            "404/",
            default_views.page_not_found,
            kwargs={"exception": Exception("Page not Found")},
        ),
        path("500/", default_views.server_error),
    ]
    if "debug_toolbar" in settings.INSTALLED_APPS:
        import debug_toolbar

        urlpatterns = [path("__debug__/", include(debug_toolbar.urls))] + urlpatterns
