from django.conf.urls.defaults import *
from django.contrib.comments.models import Comment

# Uncomment the next two lines to enable the admin:
from django.contrib import admin
admin.autodiscover()

urlpatterns = patterns('',
    # Example:
    # (r'^igem/', include('igem.foo.urls')),

    # Uncomment the admin/doc line below and add 'django.contrib.admindocs'
    # to INSTALLED_APPS to enable admin documentation:
    # (r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    (r'^admin/', include(admin.site.urls)),
    # static media
    (r'^site_media/(?P<path>.*)$', 'django.views.static.serve',{'document_root':'/igem/media'}),

    # site index
    (r'^$' , 'igem.strain_designer.views.index'),
    (r'^imptools/$' , 'igem.imptools.views.welcome'), #don 10/15/09
    (r'^imptools/start/$', 'igem.imptools.views.start'), #don 10/15/09
    (r'^imptools/start/allresults/$', 'igem.imptools.views.allresults'), #don 10/15/09
    (r'^compound/([C]\d{5})/$', 'igem.imptools.views.compound'),
    (r'^reaction/([R]\d{5})/$', 'igem.imptools.views.reaction'),
    (r'^imptools/start/allresults/detail/$', 'igem.imptools.views.detail'),
    (r'^imptools/start/allresults/detail/enzyme/$', 'igem.imptools.views.enzyme'),
    (r'^imptools/xgmml/$', 'imptools.igem.views.xgmml'),
    (r'^imptools/aboutteam/$', 'igem.imptools.views.team'),
    (r'^imptools/aboutimp/$', 'igem.imptools.views.imp'),
    (r'^imptools/demo/$', 'igem.imptools.views.demo'),
    (r'^imptools/feedback/$', 'igem.imptools.views.feedback'),
    (r'^imptools/csv/$', 'igem.imptools.views.csv'),
    (r'^imptools/bio/$', 'igem.imptools.views.prebiobrick'),
    (r'^imptools/biobrick/(?P<inputmethod>\d+)/$', 'igem.imptools.views.biobricker'),

    (r'^strain_designer/start/$', 'igem.strain_designer.views.start'),
    (r'^strain_designer/pathway/$', 'igem.strain_designer.views.imptools_results'),
    (r'^strain_designer/pathway/detail/$', 'igem.imptools.views.detail'),
    (r'^strain_designer/plasmid/detail/$', 'igem.imptools.views.detail'),
    (r'^strain_designer/pathway/detail/enzyme/$','igem.imptools.views.enzyme'),
    (r'^strain_designer/plasmid/detail/enzyme/$', 'igem.imptools.views.enzyme'),
    (r'^strain_designer/plasmid/$', 'igem.strain_designer.views.plasmid_results'),
    (r'^strain_designer/analyze/$', 'igem.strain_designer.views.analyze_results'),
    (r'^strain_designer/aboutteam/$', 'igem.strain_designer.views.team'),
    (r'^strain_designer/aboutimp/$', 'igem.strain_designer.views.about'),
    (r'^strain_designer/demo/$', 'igem.strain_designer.views.demo'),
    (r'^strain_designer/feedback/$', 'igem.strain_designer.views.feedback'),
    (r'^strain_designer/thanks/$', 'igem.strain_designer.views.feedback_thanks'),
    (r'^igem/xgmml/$', 'igem.imptools.views.xgmml'),
    (r'^igem/csv/$', 'igem.imptools.views.csv'),

    (r'^pdf_report/$', 'igem.strain_designer.views.pdf_report'),
    (r'^email_cobra/$', 'igem.strain_designer.views.email_analysis'),
    (r'^fasta_sequence/(?P<plasmid_id>\d+)/(?P<prom_id>\d+)/(?P<term_id>\d+)/(?P<rbs_id>\d+)/(?P<kegg_gene_id>[^/]+)$', 'igem.strain_designer.views.get_sequence'),
    (r'^serve_image/(?P<plasmid_id>\d+)/(?P<prom_id>\d+)/(?P<term_id>\d+)/(?P<rbs_id>\d+)/(?P<kegg_gene_id>[^/]+)$', 'igem.strain_designer.views.serve_image'),


    # ajax
    (r'^compound_name/$', 'igem.strain_designer.views.compound_name_lookup'),
    (r'^plasmid_name/$', 'igem.strain_designer.views.plasmid_name_lookup'),
    (r'^promoter_name/$', 'igem.strain_designer.views.promoter_name_lookup'),
    (r'^terminator_name/$', 'igem.strain_designer.views.terminator_name_lookup'),
    (r'^rbs_name/$', 'igem.strain_designer.views.rbs_name_lookup')

)
