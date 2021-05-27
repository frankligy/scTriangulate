"use strict";


function display_score(key,cluster) {
    // get score
    var score = document.getElementsByClassName('pass_to_js')[0].getElementsByTagName('p')[0].innerText;
    var data_from_json = JSON.parse(score);
    console.log(data_from_json);
    var display_h2_span = `${key}, ${cluster}:`;
    var display_p1 = `Reassign score is: ${data_from_json[key][0][cluster]}`;
    var display_p2 = `tf-idf score is: ${data_from_json[key][1][cluster]}`;
    var display_p3 = `SCCAF score is: ${data_from_json[key][2][cluster]}`;
    var doublet = `Average doublet score is: ${data_from_json[key][3][cluster]}`;

    var h2_obj = document.getElementById('score_information').getElementsByTagName('span');
    h2_obj[0].innerText = display_h2_span;
    var p_obj = document.getElementById('score_information').getElementsByTagName('p');
    p_obj[0].innerText = display_p1;
    p_obj[1].innerText = display_p2;
    p_obj[2].innerText = display_p3;
    var doublet_obj = document.getElementById('doublet').getElementsByTagName('h2');
    doublet_obj[0].innerText = doublet;

}

function display_plot(key,cluster) {

    var path = `./${key}_${cluster}_identity_umap.png`;
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




