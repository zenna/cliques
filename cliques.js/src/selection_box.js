cliques.SelectionBox = function (domElement) {

    // target.position is modified when panning
    this.keys =  [ 65 /*A*/, 83 /*S*/, 68 /*D*/ ];
    var box = $('<div></div>').css('position', 'fixed').css('border','1px black dashed').css('z-index','100').addClass('box').appendTo('body').hide();

    this.domElement = domElement || document;

    //internals
    var _keyPressed = false;
    var _state = this.STATE.NONE;
    var startX, startY, endX, endY;

    this.handleEvent = function ( event ) {
        if ( typeof this[ event.type ] == 'function' ) {

            this[ event.type ]( event );

        }

    };

    this.getMouseOnScreen = function( clientX, clientY ) {

        return new THREE.Vector2(
        ( clientX - this.screen.offsetLeft ) / this.radius * 0.5,
        ( clientY - this.screen.offsetTop ) / this.radius * 0.5
        );

    };

    function mousedown(event) {

        event.preventDefault();
        event.stopPropagation();

        if ( _state === this.STATE.NONE ) {

            startX = event.clientX;
            startY = event.clientY;
            box.show();
            box.css('width', 0).css('height', 0);
            box.css('left', startX);
            box.css('top', startY);

            _state = event.button;

            if ( _state === this.STATE.DRAG ) {
                console.log('start dragging');
                //_rotateStart = _rotateEnd = this.getMouseProjectionOnBall( event.clientX, event.clientY );

            }
        }
    };

    function mousemove( event ) {
        event.preventDefault();
        event.stopPropagation();

        if ( _state === this.STATE.NONE ) {

            return;

        } else if ( _state === this.STATE.DRAG ) {
            box.show();
            var mouseX = event.clientX;
            var mouseY = event.clientY;
            var w = mouseX - startX;
            var h = mouseY - startY;
            var width = Math.abs(mouseX - startX) - 5;
            var height = Math.abs(mouseY - startY) - 5;
            if (w < 0) {
                box.css('left', startX + w);
            } else {
                box.css('left', startX);
            }
            if (h < 0) {
                box.css('top', startY + h);
            } else {
                box.css('top', startY);
            }
            box.css('width', width).css('height',height);
        }

    };

    function mouseup( event ) {

        event.preventDefault();
        event.stopPropagation();

        _state = this.STATE.NONE;
        box.hide();
        console.log('stop dragging')

    };

    function bind( scope, fn ) {

        return function () {

            fn.apply( scope, arguments );

        };

    };

    // this.domElement.addEventListener( 'contextmenu', function ( event ) {
        // event.preventDefault();
    // }, false );

    this.domElement.addEventListener( 'mouseup',   bind( this, mouseup ), false );
    this.domElement.addEventListener( 'mousemove', bind( this, mousemove ), false );
    this.domElement.addEventListener( 'mousedown', bind( this, mousedown ), false );

    // window.addEventListener( 'keydown', bind( this, keydown ), false );
    // window.addEventListener( 'keyup',   bind( this, keyup ), false );

};

cliques.SelectionBox.prototype.STATE = {
    NONE : -1,
    DRAG : 0
};