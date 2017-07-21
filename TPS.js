/*
 * Matrix allocation. This function allow to create a multi-dimentional Array().
 */
function createArray( length )
{
    var arr = new Array( length || 0 );
    var i = length;
    if( arguments.length > 1 )
    {
        var args = Array.prototype.slice.call( arguments, 1 );
        while( i-- )
            arr[ length - 1 - i ] = createArray.apply( this, args );
    }
    return arr;
}
/*
 * Unit vector
 */
function unit_vector( vector )
{
    var n = vector.length;

    var norm = 0;
    for( var i = 0; i < n; i++ )
        norm += vector[ i ] * vector[ i ];
    norm = Math.sqrt( norm );

    var r = new Array();
    for( var i = 0; i < n; i++ )
        r[ i ] = vector[ i ] / norm;

    return r;
}

/*
 * Functions defined in the publication of Bookstein (1989). The original
 * function U is not realy used here because of performance issues. Since the U
 * function is used for the K matrix, the sqrt() present in the distance
 * function is simplified with the ^2. The euclidean distance being neved used
 * elsewhere, it's not present in the TPS library.
 */
function U( r )
{
    return U2( r * r );
}

function U2( r )
{
    if( r == 0 )
        return 0;
    else
        return r * Math.log( r );
}

/*
 * The generate function allow to calculate the TPS parameters of the
 * distortion. The parameters 'src' (source of the transoformation) and 'dst'
 * (destination points) are mendatory.
 */
function TPS_generate( src, dst )
{
    var n = src.length;

    // K Matrix
    var K = createArray( n, n );

    for( var i = 0; i < n; i++ )
    {
        for( var j = 0; j < n; j++ )
        {
            var dx = src[ j ][ 0 ] - src[ i ][ 0 ];
            var dy = src[ j ][ 1 ] - src[ i ][ 1 ];
            K[ i ][ j ] = U2( dx * dx + dy * dy );
        }
    }

    // L Matrix
    var L = createArray( n + 3, n + 3 );

    for( var i = 0; i < n; i++ )
    {
        for( var j = 0; j < n; j++ )
        {
            L[ i ][ j ] = K[ i ][ j ]
        }
    }

    for( var i = 0; i < n; i++ )
    {
        L[ i ][ n ] = 1;
        L[ i ][ n + 1 ] = src[ i ][ 0 ];
        L[ i ][ n + 2 ] = src[ i ][ 1 ];

        L[ n ][ i ] = 1;
        L[ n + 1 ][ i ] = src[ i ][ 0 ];
        L[ n + 2 ][ i ] = src[ i ][ 1 ];
    }

    for( var i = 0; i < 3; i++ )
    {
        for( var j = 0; j < 3; j++ )
        {
            L[ n + i ][ n + j ] = 0;
        }
    }

    // V Matrix
    var V = createArray( n + 3, 2 );
    for( var i = 0; i < n; i++ )
    {
        V[ i ][ 0 ] = dst[ i ][ 0 ];
        V[ i ][ 1 ] = dst[ i ][ 1 ];
    }

    for( var i = 0; i < 3; i++ )
    {
        V[ i + n ][ 0 ] = 0;
        V[ i + n ][ 1 ] = 0;
    }

    // Solve
    var L = new Matrix( L );
    var V = new Matrix( V );

    var Wa = L.solve( V ).toArray();

    var W = Wa.slice( 0, n );
    var a = Wa.slice( n, n + 3 );

    // Return
    var g = {
    src : src,
    dst : dst,
    linear : a,
    weights : W
    }

    return g;
}

/*
 * Projection function. This function take as parameter the 'g' TPS parameters
 * (see the 'generate' function), and the 'x' and 'y' coordinates of the point
 * to project.
 */
function TPS_project( g, x, y )
{
    // vars
    var n = g[ 'src' ].length

    var xy = new Array();
    xy[ 0 ] = 1;
    xy[ 1 ] = x;
    xy[ 2 ] = y;

    // Linear part -- dot( hstack( ( 1, XY ) ), linear )
    var p = new Array();
    for( var i = 0; i < 2; i++ )
    {
        var tmp = 0;
        for( var j = 0; j < 3; j++ )
        {
            tmp += xy[ j ] * g[ 'linear' ][ j ][ i ];
        }

        p[ i ] = tmp;
    }

    // Distance with the src vector
    var dist = new Array();
    for( var i = 0; i < n; i++ )
    {
        var tmp = 0;
        // sum( ( src - XY ) ** 2
        for( var j = 0; j < 2; j++ )
        {
            tmp += Math.pow( ( g[ 'src' ][ i ][ j ] - xy[ j + 1 ] ), 2 );
        }

        // U2
        dist[ i ] = U2( tmp );
    }

    // Non-linear part -- dot( ... , W )
    for( var i = 0; i < 2; i++ )
    {
        for( var j = 0; j < n; j++ )
        {
            p[ i ] += dist[ j ] * g[ 'weights' ][ j ][ i ];
        }
    }

    // Return
    return p;
}
