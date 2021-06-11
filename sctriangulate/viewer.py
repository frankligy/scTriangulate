from yattag import Doc
import json

# for individual viewer
def html_banner():
    doc,tag,text,line = Doc().ttl()
    with tag('div'):
        doc.attr(id='banner')
        line('h2','scTriangulate Viewer')
    return doc.getvalue()


def html_left_nav(key_cluster_dict):
    doc,tag,text,line = Doc().ttl()
    with tag('div'):
        doc.attr(id='left_nav')
        for key,value in key_cluster_dict.items():
            line('h3',key)
            with tag('ul'):
                doc.attr(id='{}'.format(key))
                for cluster in value:
                    line('li',cluster,onclick='display_score(this.parentElement.id,this.textContent);display_plot(this.parentElement.id,this.textContent)',
                    onmouseover='change_color(this)',onmouseout='change_color_back(this)')
    return doc.getvalue()

def html_right_show(key_cluster_data,total_metrics):
    doc,tag,text,line = Doc().ttl()
    with tag('div'):
        doc.attr(id='right_show')
        with tag('div'):
            doc.attr(id='score_information')
            with tag('h2'):
                text('Scoring information:')
                with tag('span'):
                    text(' ')
            with tag('table'):
                for metric in total_metrics:
                    with tag('tr'):
                        line('td',metric)
                        line('td','',id='cluster_to_{}'.format(metric))

        with tag('div'):
            doc.attr(id='identity')
            line('h2','Cluster location')
            doc.stag('img',src='./init.png',width='60%',height='60%',alt='Choose key and cluster',klass='img_identity')                
        
        with tag('div'):
            doc.attr(id='doublet')
            line('h2','doublet distribution')
            doc.stag('img',src='./umap_sctriangulate_doublet_scores.png',width='60%',height='60%')

        with tag('div'):
            doc.attr(id='enrichment')
            line('h2','Enrichment plot')
            doc.stag('img',src='./init.png',alt='Choose key and cluster',klass='img_enrichment',width='60%',height='60%')
        with tag('div'):
            doc.attr(klass='marker_umap')
            line('h2','Marker gene umap')
            doc.stag('img',src='./init.png',alt='Choose key and cluster',klass='img_marker',width='100%')
        with tag('div'):
            doc.attr(klass='exclusive_umap')
            line('h2','Exclusive gene umap')
            doc.stag('img',src='./init.png',alt='Choose key and cluster',klass='img_exclusive',width='100%')
        with tag('div'):
            doc.attr(klass='inspection_div')
            line('a','Go to Inspection',href='./inspection.html')
        with tag('div'):
            doc.attr('hidden',klass='pass_to_js')
            with tag('p'):
                text(json.dumps(key_cluster_data))
    return doc.getvalue()


def to_html(key_cluster_dict,key_cluster_data,total_metrics):
    doc,tag,text,line = Doc().ttl()
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            with tag('title'):
                text('scTriangulate_individual')
            with tag('link'):
                doc.attr(type='text/css',rel='stylesheet',href='./viewer.css')
        with tag('body'):
            doc.asis(html_banner())
            doc.asis(html_left_nav(key_cluster_dict))
            doc.asis(html_right_show(key_cluster_data,total_metrics))
        with tag('script'):
            doc.attr(src='./viewer.js')
    return doc.getvalue()



# for inspection viewer
def header():
    doc,tag,text,line = Doc().ttl()
    line('h1','ScTriangulate Inspection',id='header')
    return doc.getvalue()

def left_nav(key_cluster_dict,reference):
    doc,tag,text,line = Doc().ttl()
    data = key_cluster_dict[reference] # a list
    with tag('div'):
        doc.attr(id='left_nav')
        with tag('h3'):
            text('{}'.format(reference))
        with tag('ul'):
            doc.attr(id=reference)
            for cluster in data:
                line('li',cluster,onclick='display(this.parentElement.id,this.textContent)',onmouseover='change_color(this)',onmouseout='change_color_back(this)')
    return doc.getvalue()

def right_show():
    doc,tag,text,line = Doc().ttl()
    with tag('div'):
        doc.attr(id='right_show')
        line('h2','Cluster',id='identity')
        with tag('div'):
            doc.attr(id='two_umap_left')
            with tag('div'):
                doc.attr(klass='umap_div')
                line('h2','UMAP view')
                doc.stag('img', id='umap', src='./init.png', alt='please check cluster on the left')
        with tag('div'):
            doc.attr(klass='heatmap_div')
            line('h2','heatmap view')
            doc.stag('img', id='heatmap', src='./init.png', alt='please check cluster on the left')
        with tag('div'):
            doc.attr(klass='viewer_div')
            line('a','Go to Viewer',href='./viewer.html')
    return doc.getvalue()


def inspection_html(key_cluster_dict,reference):
    doc,tag,text,line = Doc().ttl()
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            doc.stag('meta', charset='UTF-8')
            doc.stag('meta', name='viewport', content='width=device-width, initial-scale=1.0')
            with tag('title'):
                text('scTriangulate_inspection')
            with tag('link'):
                doc.attr(type='text/css',rel='stylesheet',href='./inspection.css')
        with tag('body'):
            doc.asis(header())
            doc.asis(left_nav(key_cluster_dict,reference))
            doc.asis(right_show())
        with tag('script'):
            doc.attr(src='./inspection.js')
    return doc.getvalue()






