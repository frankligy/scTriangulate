"use strict";


function display_score(key,cluster) {
    var data_from_json = JSON.parse(score);
    var display_h2 = `Score information of ${key}, ${cluster}:`;
    var display_p1 = `Reassign score is: ${data_from_json.key[0].cluster}`;
    var display_p2 = `tf-idf score is: ${data_from_json.key[1].cluster}`;
    var display_p3 = `SCCAF score is: ${data_from_json.key[2].cluster}`;

    var h2_obj = document.getElementsByClassName('score_h2');
    h2_obj[0].innerText = display_h2;
    var p1_obj = document.getElementsByClassName('score_p1');
    p1_obj[0].innerText = display_p1;
    var p2_obj = document.getElementsByClassName('score_p2');
    p2_obj[0].innerText = display_p2;
    var p3_obj = document.getElementsByClassName('score_p3');
    p2_obj[0].innerText = display_p3;

}

function display_plot(key,cluster) {
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




