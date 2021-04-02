#!/data/salomonis2/LabFiles/Frank-Li/citeseq/scanpy_new_env/bin/python3.6

from yattag import Doc
import json


def html_banner():
    doc,tag,text,line = Doc().ttl()
    with tag('div'):
        doc.attr(klass='banner')
        with tag('h2'):
            text('scTriangulate Viewer')
    return doc.getvalue()


def html_left_nav(key_cluster_dict):
    doc,tag,text,line = Doc().ttl()
    with tag('div'):
        doc.attr(klass='left_nav')
        for key,value in key_cluster_dict.items():
            with tag('h3'):
                text(key)
            with tag('ul'):
                doc.attr(klass='{}'.format(key))
                for cluster in value:
                    with tag('li'):
                        doc.attr(onclick='display_score(this.parentElement.className,this.textContent);display_plot(this.parentElement.className,this.textContent)')
                        text(cluster)
    return doc.getvalue()

def html_right_show(key_cluster_data):
    doc,tag,text,line = Doc().ttl()
    with tag('div'):
        doc.attr(klass='right_show')
        with tag('div'):
            doc.attr(klass='score_information')
            with tag('h2'):
                doc.attr(klass='score_h2')
                text('placeholder')
            with tag('p'):
                doc.attr(klass='score_p1')
                text('placeholder')
            with tag('p'):
                doc.attr(klass='score_p2')
                text('placeholder')   
            with tag('p'):
                doc.attr(klass='score_p3')
                text('placeholder')  
        with tag('div'):
            doc.attr(klass='doublet')
            with tag('h2'):
                doc.attr(klass='doublet_h2')
                text('placeholder')
            doc.stag('img',src='./doublet.png',width='60%')

        with tag('div'):
            doc.attr(klass='enrichment')
            with tag('h2'):
                text('enrichment plot')
            doc.stag('img',src='./init.png',alt='Choose key and cluster',klass='img_enrichment',width='60%')
        with tag('div'):
            doc.attr(klass='marker_umap')
            with tag('h2'):
                text('marker gene umap')
            doc.stag('img',src='./init.png',alt='Choose key and cluster',klass='img_marker',width='90%')
        with tag('div'):
            doc.attr(klass='exclusive_umap')
            with tag('h2'):
                text('exclusive gene umap')
            doc.stag('img',src='./init.png',alt='Choose key and cluster',klass='img_exclusive',width='90%')

        with tag('div'):
            doc.attr(klass='pass_to_js')
            with tag('p'):
                text(json.dumps(key_cluster_data))
    return doc.getvalue()


def to_html(key_cluster_dict,key_cluster_data):
    doc,tag,text,line = Doc().ttl()
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            with tag('title'):
                text('scTriangulate')
            with tag('link'):
                doc.attr(rel='stylesheet',href='./viewer.css')
        with tag('body'):
            doc.asis(html_banner())
            with tag('div'):
                doc.attr(klass='main_body')
                doc.asis(html_left_nav(key_cluster_dict))
                doc.asis(html_right_show(key_cluster_data))
        with tag('script'):
            doc.attr(src='../score.json')
            doc.attr(src='./viewer.js')
    return doc.getvalue()







