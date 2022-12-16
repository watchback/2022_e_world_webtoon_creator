/*여기부터 SVG 드래그 관련 Js코드임*/
function makeDraggable(evt) {
    var selectedElement = false;
    var svg = evt.target;
    var selectedElement, offset;
    var mx, my;

    svg.addEventListener('mousedown', startDrag);
    svg.addEventListener('mousemove', drag);
    svg.addEventListener('mouseup', endDrag);
    svg.addEventListener('mouseleave', endDrag);


    function getMousePosition(evt) {
        var CTM = svg.getScreenCTM();
        return {
            cx: (evt.clientX - CTM.e) / CTM.a,
            cy: (evt.clientY - CTM.f) / CTM.d
        };
    }
    function startDrag(evt) {
        if (evt.target.classList.contains('draggable')) {
            selectedElement = evt.target;
            offset = getMousePosition(evt);
            offset.cx -= parseFloat(selectedElement.getAttributeNS(null, "cx"));
            offset.cy -= parseFloat(selectedElement.getAttributeNS(null, "cy"));
        }
        if (evt.target.classList.contains('static')) {
            selectedElement = evt.target;
            offset = getMousePosition(evt);
            offset.x -= parseFloat(selectedElement.getAttributeNS(null, "cx"));
            offset.y -= parseFloat(selectedElement.getAttributeNS(null, "cy"));
        }
    }
    function drag(evt) {
        if (selectedElement) {
            evt.preventDefault();
            var coord = getMousePosition(evt);
            selectedElement.setAttributeNS(null, "cx", coord.cx - offset.cx);
            selectedElement.setAttributeNS(null, "cy", coord.cy - offset.cy);
        }
    }
    function endDrag(evt) {
        selectedElement = null;
    }
}