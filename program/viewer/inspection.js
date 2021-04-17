"use strict";

function display(cluster) {
    var h2_div = document.getElementById('identity');
    h2_div.innerText = cluster;

    var path = `./UMAP_${cluster}.pdf`;
    var umap_div = document.getElementById('umap');
    umap_div.src = path;

    var path = `./DE_heatmap_${cluster}.pdf`;
    var heatmap_div = document.getElementById('heatmap');
    heatmap_div.src = path;
}