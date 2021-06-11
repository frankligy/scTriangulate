"use strict";


function display_score(key,cluster) {
    // get score
    var score = document.getElementsByClassName('pass_to_js')[0].getElementsByTagName('p')[0].innerText;
    var data_from_json = JSON.parse(score);
    var display_h2_span = `${key}, ${cluster}:`;
    var h2_obj = document.getElementById('score_information').getElementsByTagName('span');
    h2_obj[0].innerText = display_h2_span;
    var scores = Object.keys(data_from_json[key]);
    var scores_len = scores.length;
    console.log(scores)
    for (let i=0; i<scores_len; i++) {
        let display = `${data_from_json[key][scores[i]][cluster]}`;
        let display_obj = document.getElementById(scores[i]);
        display_obj.innerText = display;
    }

}

function display_plot(key,cluster) {

    var path = `./${key}_${cluster}_location_umap.png`;
    var img_identity = document.getElementsByClassName('img_identity');
    img_identity[0].src = path;

    var path = `./${key}_${cluster}_enrichment.png`;
    var img_enrichment = document.getElementsByClassName('img_enrichment');
    img_enrichment[0].src = path;

    var path = `./${key}_${cluster}_marker_umap.png`;
    var img_marker = document.getElementsByClassName('img_marker');
    img_marker[0].src = path;

    var path = `./${key}_${cluster}_exclusive_umap.png`;
    var img_marker = document.getElementsByClassName('img_exclusive');
    img_marker[0].src = path;
}

function change_color(object) {
    object.style.color = 'blue';
}

function change_color_back(object) {
    object.style.color = 'black';
}




