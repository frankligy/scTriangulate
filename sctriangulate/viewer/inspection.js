"use strict";

function display(key,cluster) {
    var h2_div = document.getElementById('identity');
    h2_div.innerText = cluster;

    var path = `./${cluster}_heterogeneity_umap.png`;
    var umap_div = document.getElementById('umap');
    umap_div.src = path;

    var path = `./${cluster}_heterogeneity_heatmap.pdf`;
    var heatmap_div = document.getElementById('heatmap');
    heatmap_div.src = path;
}

function change_color(object) {
    object.style.color = 'blue';
}

function change_color_back(object) {
    object.style.color = 'black';
}
