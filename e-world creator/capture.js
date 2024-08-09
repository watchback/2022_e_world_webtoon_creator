function click(){
    html2canvas(document.querySelector('.canvas-container'), {}).then(function (canvas) {
        saveAs(canvas.toDataURL(), 'webtoonsaved.png');
    });    
}

function saveAs(uri, filename) {
	var link = document.createElement('a');
	if (typeof link.download === 'string') {
		link.href = uri;
		link.download = filename;
		document.body.appendChild(link);
		link.click();
		document.body.removeChild(link);
	} else {
		window.open(uri);
	}
}