/**
 * Code generating polygons set for rendering beachballs in 2D space.
 * It's a TypeScript version of patch_meca.m file found in Mirone source code
 * (http://w3.ualg.pt/~jluis/mirone/) - a MATLAB-based framework distrubuted
 * under LGPL license.
 *
 * Comments within this function borrowed from the .m file.
 *
 * Generated polygons coordinates are located within following axes range:
 * x axis -> [-1,1]
 * y axis -> [-1,1]
 * So you need to translate and scale coordinates according to the need of your
 * rendering code.
 *
 * NOTE: code layout is not optimal to perserve mapping to MATLAB script in
 * Mirone
 *
 */

import * as bbmath from './bbmath'

var sylvester = require('sylvester')

class SDR {
    str: number
    dip: number
    rake: number
}

var np = 45
var epsilon = 1E-5

var sign = (x: number) => x > 0 ? 1 : x < 0 ? -1 : 0

var range = function (start: number, end: number, step: number) {
    var range: number[] = [];
    while (step > 0 ? end >= start : end <= start) {
        range.push(start);
        start += step;
    }
    return range;
};

/**
 * Put an angle between 0 and 360 degrees
 * @param {type} ang
 * @returns {undefined}
 */
var zero_360 = (ang: number) => {
    let str: number
    if (ang >= 360.0) {
        str = ang - 360.0
    } else if (ang < 0.0) {
        str = ang + 360.0
    } else {
        str = ang;
    }
    return str;
}

/**
 * Compute arctg in degrees, between -180 et +180.
 * @param y
 * @param x
 * @returns
 */
var datan2 = (y: number, x: number) => {
    var arctg: number
    if (Math.abs(x) < bbmath.EPS) {
        if (Math.abs(y) < bbmath.EPS) {
            arctg = 0.0
        }
        else if (y < 0) {
            arctg = -90
        }
        else {
            arctg = 90
        }
    } else if (x < 0) {
        if (y < 0) {
            arctg = Math.atan(y / x) * bbmath.R2D - 180
        } else {
            arctg = Math.atan(y / x) * bbmath.R2D + 180
        }
    } else {
        arctg = Math.atan(y / x) * bbmath.R2D
    }
    return arctg
}

/**
 * Compute the strike of the decond nodal plane when are given
 * strike, dip and rake for the first nodal plane with AKI & RICHARD's
 * convention. Angles are in degrees.
 * @param NP1
 * @param NP2
 */
var compute_strike = (NP1: SDR, NP2: SDR) => {
    var sr = Math.sin(NP1.rake * bbmath.D2R)
    var cr = Math.cos(NP1.rake * bbmath.D2R)
    var ss = Math.sin(NP1.str * bbmath.D2R)
    var cs = Math.cos(NP1.str * bbmath.D2R)
    var cd1 = Math.cos(NP1.dip * bbmath.D2R)
    var am: number, temp: number, sp2: number, cp2: number, str2: number
    if (NP1.rake === 0) {
        am = 1
    } else {
        am = NP1.rake / Math.abs(NP1.rake)
    }

    if (cd1 < bbmath.EPS && Math.abs(cr) < bbmath.EPS) {
        NP2.str = NP1.str + 180
    } else {
        temp = cr * cs
        temp = temp + sr * ss * cd1
        sp2 = -am * temp
        temp = ss * cr
        temp = temp - sr * cs * cd1
        cp2 = am * temp
        str2 = datan2(sp2, cp2)
        str2 = zero_360(str2)
        NP2.str = str2
    }
}

/**
 * Computes and sets second nodal plane dip when strike,
 * dip and rake for the first nodal plane specified with AKI & RICHARD's
 * convention. Angles are in degrees.
 *
 * @param {SDR} NP1 - first nodal plane
 * @param {SDR} NP2 - second nodal plane
 */
var compute_dip = (NP1: SDR, NP2: SDR) => {
    var am: number = NP1.rake === 0 ? 1 : (NP1.rake / Math.abs(NP1.rake))
    NP2.dip = Math.acos(am * Math.sin(NP1.rake * bbmath.D2R) * Math.sin(NP1.dip * bbmath.D2R)) / (bbmath.D2R)
}

/**
 * Compute rake in the second nodal plane when strike ,dip
 * and rake are given for the first nodal plane with AKI &
 * RICHARD's convention. Angles are in degrees.
 *
 * @param {SDR} NP1
 * @param {SDR} NP2
 */
var compute_rake = (NP1: SDR, NP2: SDR) => {
    var str2 = NP2.str
    var dip2 = NP2.dip

    var am: number = NP1.rake === 0 ? 1 : (NP1.rake / Math.abs(NP1.rake))
    var sd = Math.sin(NP1.dip * bbmath.D2R)
    var cd = Math.cos(NP1.dip * bbmath.D2R)
    var ss = Math.sin((NP1.str - str2) * bbmath.D2R)
    var cs = Math.cos((NP1.str - str2) * bbmath.D2R)

    var sinrake2 = Math.abs(dip2 - 90.0) < bbmath.EPS ? (am * cd) : (-am * sd * cs / cd)

    NP2.rake = datan2(sinrake2, -am * sd * ss)
}

/**
 * Computes second nodal plane
 *
 * @param {SDR}
 * @return {SDR} - second nodal plane
 */
var define_second_plane = (NP1: SDR): SDR => {
    var NP2 = {
        str: 0,
        dip: 0,
        rake: 0
    }
    compute_strike(NP1, NP2)
    compute_dip(NP1, NP2)
    compute_rake(NP1, NP2)
    return NP2
}

/**
 * Compute null axis strike when strike and dip are given
 * for each nodal plane. Angles are in degrees.
 *
 * @param {SDR} NP1
 * @param {SDR} NP2
 * @returns {number} null axis strike angle
 */
var null_axis_strike = (NP1: SDR, NP2: SDR): number => {
    var sd1: number, sd2: number, cd1: number, cd2: number, ss1: number, ss2: number, cs1: number, cs2: number, cosphn: number, sinphn: number, phn: number

    sd1 = Math.sin(NP1.dip * bbmath.D2R)
    cd1 = Math.cos(NP1.dip * bbmath.D2R)
    sd2 = Math.sin(NP2.dip * bbmath.D2R)
    cd2 = Math.cos(NP2.dip * bbmath.D2R)
    ss1 = Math.sin(NP1.str * bbmath.D2R)
    cs1 = Math.cos(NP1.str * bbmath.D2R)
    ss2 = Math.sin(NP2.str * bbmath.D2R)
    cs2 = Math.cos(NP2.str * bbmath.D2R)

    cosphn = sd1 * cs1 * cd2 - sd2 * cs2 * cd1
    sinphn = sd1 * ss1 * cd2 - sd2 * ss2 * cd1

    if (Math.sin((NP1.str - NP2.str) * bbmath.D2R) < 0) {
        cosphn = -cosphn
        sinphn = -sinphn
    }

    phn = datan2(sinphn, cosphn)
    phn = phn < 0 ? phn + 360 : phn
    return phn
}

/**
 * Compute null axis dip when strike and dip are given
 * for each nodal plane. Angles are in degrees.
 *
 * @param {SDR} NP1
 * @param {SDR} NP2
 * @returns {number} null axis dip angle
 */
var null_axis_dip = (NP1: SDR, NP2: SDR) => {
    var dip = Math.asin(Math.sin(NP1.dip * bbmath.D2R) * Math.sin(NP2.dip * bbmath.D2R) * Math.sin((NP1.str - NP2.str) * bbmath.D2R)) / bbmath.D2R
    dip = dip < 0 ? -dip : dip
    return dip
}

/**
 * From FORTRAN routines of Anne Deschamps :
 * compute azimuth and plungement of P-T axis
 * from nodal plane strikes, dips and rakes.
 *
 * @param {SDR} NP1
 * @param {SDR} NP2
 */
var dc_to_axe = (NP1: SDR, NP2: SDR) => {
    var P: SDR = new SDR()
    var T: SDR = new SDR()
    var N: SDR = new SDR()
    var im: number, cd1: number, cd2: number, sd1: number, sd2: number, cp1: number, cp2: number, sp1: number, sp2: number,
        amz: number, amx: number, amy: number, dx: number, px: number, dy: number, py: number
    var pure_strike_slip = 0

    if (Math.abs(Math.sin(NP1.rake * bbmath.D2R)) > bbmath.EPS) {
        im = NP1.rake / Math.abs(NP1.rake)
    } else if (Math.abs(Math.sin(NP2.rake * bbmath.D2R)) > bbmath.EPS) {
        im = NP2.rake / Math.abs(NP2.rake)
    } else {
        pure_strike_slip = 1
    }

    if (pure_strike_slip) {
        if (Math.cos(NP1.rake * bbmath.D2R) < 0.0) {
            P.str = zero_360(NP1.str + 45)
            T.str = zero_360(NP1.str - 45)
        } else {
            P.str = zero_360(NP1.str - 45)
            T.str = zero_360(NP1.str + 45)
        }
        P.dip = 0
        T.dip = 0
    } else {
        cd1 = Math.cos(NP1.dip * bbmath.D2R) * bbmath.MSQRT2
        sd1 = Math.sin(NP1.dip * bbmath.D2R) * bbmath.MSQRT2
        cd2 = Math.cos(NP2.dip * bbmath.D2R) * bbmath.MSQRT2
        sd2 = Math.sin(NP2.dip * bbmath.D2R) * bbmath.MSQRT2
        cp1 = -Math.cos(NP1.str * bbmath.D2R) * sd1
        sp1 = Math.sin(NP1.str * bbmath.D2R) * sd1
        cp2 = -Math.cos(NP2.str * bbmath.D2R) * sd2
        sp2 = Math.sin(NP2.str * bbmath.D2R) * sd2

        amz = -(cd1 + cd2)
        amx = -(sp1 + sp2)
        amy = cp1 + cp2
        dx = Math.atan2(Math.sqrt(amx * amx + amy * amy), amz) - bbmath.MPI2
        px = Math.atan2(amy, -amx)
        px = px < 0 ? px + bbmath.TWO_PI : px

        amz = cd1 - cd2
        amx = sp1 - sp2
        amy = -cp1 + cp2
        dy = Math.atan2(Math.sqrt(amx * amx + amy * amy), -Math.abs(amz)) - bbmath.MPI2
        py = Math.atan2(amy, -amx)
        py = amz > 0 ? py - Math.PI : py
        py = py < 0 ? py + bbmath.TWO_PI : py

        if (im === 1) {
            P.dip = dy
            P.str = py
            T.dip = dx
            T.str = px
        } else {
            P.dip = dx
            P.str = px
            T.dip = dy
            T.str = py
        }
    }

    T.str = T.str / bbmath.D2R
    T.dip = T.dip / bbmath.D2R
    P.str = P.str / bbmath.D2R
    P.dip = P.dip / bbmath.D2R

    N.dip = null_axis_dip(NP1, NP2)
    N.str = null_axis_strike(NP1, NP2)

    return {
        P: P,
        T: T,
        N: N
    }
}

var doublecouple = function (str1: number, dip1: number, rake1: number, str2?: number, dip2?: number, rake2?: number) {
    var NP1 = {
        str: 0,
        dip: 0,
        rake: 0
    }
    var NP2 = {
        str: 0,
        dip: 0,
        rake: 0
    }
    var comp: any, dilat: any

    dip1 = Math.abs(dip1 - 90) < epsilon ? 90 - epsilon : dip1

    rake1 = rake1 > 180 ? rake1 - 360 : rake1

    if (Math.abs(rake1) - 90 < epsilon || Math.abs(rake1) - 180 < epsilon) {
        rake1 = rake1 - sign(rake1) * epsilon
    }

    if (arguments.length === 6) {
        dip2 = Math.abs(dip2 - 90) < epsilon ? dip2 = 90 - epsilon : dip2
        rake2 = rake2 > 180 ? rake2 - 360 : rake2
        if (Math.abs(rake2) - 90 < epsilon || Math.abs(rake1) - 180 < epsilon) {
            rake2 = rake2 - sign(rake2) * epsilon
        }
    }

    NP1.str = str1
    NP1.dip = dip1
    NP1.rake = rake1

    if (arguments.length === 3) {
        NP2 = define_second_plane(NP1)
    } else {
        NP2.str = str2
        NP2.dip = dip2
        NP2.rake = rake2
    }

    var PTN = dc_to_axe(NP1, NP2)
    var N = PTN.N

    // First nodal plane
    var rot_s1 = sylvester.Matrix.create([
        [Math.cos(NP1.str * bbmath.D2R), -Math.sin(NP1.str * bbmath.D2R)],
        [Math.sin(NP1.str * bbmath.D2R), Math.cos(NP1.str * bbmath.D2R)]
    ])

    var ang1 = (N.str - NP1.str) * bbmath.D2R

    ang1 = ang1 < 0 ? ang1 + bbmath.TWO_PI : ang1
    ang1 = ang1 > Math.PI ? ang1 - Math.PI : ang1

    var teta1 = range(0, ang1, Math.PI / np)
    teta1.push(ang1)

    var teta2 = range(ang1, Math.PI, Math.PI / np)
    teta2.push(Math.PI)

    var s_teta1 = teta1.map(item => Math.sin(item))
    var s_teta2 = teta2.map(item => Math.sin(item))

    var radip = s_teta1.map(item => Math.atan(Math.tan(NP1.dip * bbmath.D2R) * item))

    var rproj = radip.map(item => bbmath.MSQRT2 * Math.sin((bbmath.MPI2 - item) / 2))

    var x1_P1 = rproj.map((item, idx) => item * s_teta1[idx])

    var y1_P1 = rproj.map((item, idx) => item * Math.cos(teta1[idx]))

    // x2
    radip = s_teta2.map(item => Math.atan(Math.tan(NP1.dip * bbmath.D2R) * item))

    rproj = radip.map(item => bbmath.MSQRT2 * Math.sin((bbmath.MPI2 - item) / 2))

    var x2_P1 = rproj.map((item, idx) => item * s_teta2[idx])

    var y2_P1 = rproj.map((item, idx) => item * Math.cos(teta2[idx]))

    // Rotate to P1_strike
    var x1trans = (sylvester.Matrix.create(x1_P1)).transpose()
    var y1trans = (sylvester.Matrix.create(y1_P1)).transpose()
    var transposedMatrix: number[][] = []
    x1trans.elements[0].forEach((item: number, idx:number) => transposedMatrix.push([item, y1trans.elements[0][idx]]))

    var x1y1 = sylvester.Matrix.create(transposedMatrix)
    var plan1_a = x1y1.x(rot_s1)


    var x2trans = (sylvester.Matrix.create(x2_P1)).transpose()
    var y2trans = (sylvester.Matrix.create(y2_P1)).transpose()
    transposedMatrix = []
    x2trans.elements[0].forEach((item: number, idx: number) => transposedMatrix.push([item, y2trans.elements[0][idx]]))

    var x2y2 = sylvester.Matrix.create(transposedMatrix)
    var plan1_b = x2y2.x(rot_s1)

    // second nodal plane
    var rot_s2 = sylvester.Matrix.create([
        [Math.cos(NP2.str * bbmath.D2R), -Math.sin(NP2.str * bbmath.D2R)],
        [Math.sin(NP2.str * bbmath.D2R), Math.cos(NP2.str * bbmath.D2R)]
    ])
    var ang2 = (NP2.str - N.str) * bbmath.D2R
    ang2 = ang2 < 0 ? ang2 + bbmath.TWO_PI : ang2
    ang2 = ang2 > Math.PI ? ang2 - Math.PI : ang2

    teta1 = range(0, Math.PI - ang2, Math.PI / np)
    teta1.push(Math.PI - ang2)
    teta2 = range(Math.PI - ang2, Math.PI, Math.PI / np)
    teta2.push(Math.PI)

    s_teta1 = teta1.map(item => Math.sin(item))
    s_teta2 = teta2.map(item => Math.sin(item))

    radip = s_teta1.map(item => Math.atan(Math.tan(NP2.dip * bbmath.D2R) * item))

    rproj = radip.map(item => bbmath.MSQRT2 * Math.sin((bbmath.MPI2 - item) / 2))

    var x1_P2 = rproj.map((item, idx) => item * s_teta1[idx])

    var y1_P2 = rproj.map((item, idx) => item * Math.cos(teta1[idx]))

    // x2
    radip = s_teta2.map((item) => Math.atan(Math.tan(NP2.dip * bbmath.D2R) * item))

    rproj = radip.map((item) => bbmath.MSQRT2 * Math.sin((bbmath.MPI2 - item) / 2))

    var x2_P2 = rproj.map((item, idx) => item * s_teta2[idx])

    var y2_P2 = rproj.map((item, idx) => item * Math.cos(teta2[idx]))

    // Rotate to P1_strike
    x1trans = (sylvester.Matrix.create(x1_P2)).transpose()
    y1trans = (sylvester.Matrix.create(y1_P2)).transpose()
    transposedMatrix = []
    x1trans.elements[0].forEach((item: number, idx: number) => transposedMatrix.push([item, y1trans.elements[0][idx]]))

    var x1y2 = sylvester.Matrix.create(transposedMatrix)
    var plan2_a = x1y2.x(rot_s2)

    x2trans = (sylvester.Matrix.create(x2_P2)).transpose()
    y2trans = (sylvester.Matrix.create(y2_P2)).transpose()
    transposedMatrix = []
    x2trans.elements[0].forEach((item: number, idx: number) => transposedMatrix.push([item, y2trans.elements[0][idx]]))

    x2y2 = sylvester.Matrix.create(transposedMatrix)
    var plan2_b = x2y2.x(rot_s2)

    //Compute the arc-circles that close compressive/dilatational parts
    var a1: number, a2: number, aa: number, ang_g: number, ang_l: number
    if (NP1.rake >= 0) {
        if (NP2.str < NP1.str) {
            a1 = NP2.str * bbmath.D2R + Math.PI
            a2 = NP1.str * bbmath.D2R
            ang_g = Math.max(a1, a2)
            ang_l = Math.min(a1, a2)
        } else {
            a1 = NP2.str * bbmath.D2R - Math.PI
            a2 = NP1.str * bbmath.D2R
            ang_g = Math.max(a1, a2)
            ang_l = Math.min(a1, a2)
        }
    } else {
        ang_g = NP1.str * bbmath.D2R
        ang_l = NP2.str * bbmath.D2R
        ang_g = ang_l > ang_g ? ang_g + bbmath.TWO_PI : ang_g
    }

    var tetaRange = range(ang_l, ang_g, Math.PI / np)
    tetaRange.push(ang_g)
    var teta = tetaRange.map((item) => Math.PI / 2 - item)
    teta = teta.reverse()

    var arc1_x = teta.map(item => Math.cos(item))
    var arc1_y = teta.map(item => Math.sin(item))

    teta = teta.reverse()
    var arc3_x = teta.map(item => -Math.cos(item))
    var arc3_y = teta.map(item => -Math.sin(item))

    if (NP1.rake >= 0) {
        if (NP2.str > NP1.str) {
            ang_g = NP1.str * bbmath.D2R + Math.PI
            ang_l = NP2.str * bbmath.D2R - Math.PI
        }
        else {
            ang_g = NP1.str * bbmath.D2R + Math.PI
            ang_l = NP2.str * bbmath.D2R + Math.PI
        }

    } else {
        ang_l = NP1.str * bbmath.D2R
        ang_g = NP2.str * bbmath.D2R + Math.PI
        if (ang_g - ang_l > bbmath.TWO_PI) {
            ang_g = ang_g - bbmath.TWO_PI
        }
        if (ang_l > ang_g) {
            var tmp: number
            ang_l = NP2.str * bbmath.D2R
            ang_g = NP1.str * bbmath.D2R - Math.PI
            tmp = plan1_a
            plan1_a = plan2_a
            plan2_a = tmp
            tmp = plan1_b
            plan1_b = plan2_b
            plan2_b = tmp
        }
    }

    tetaRange = range(ang_l, ang_g, Math.PI / np)
    tetaRange.push(ang_g)
    teta = tetaRange.map(item => Math.PI / 2 - item)

    var arc2_x = teta.map(item => Math.cos(item))
    var arc2_y = teta.map(item => Math.sin(item))

    teta = teta.reverse()
    var arc4_x = teta.map(item => -Math.cos(item))
    var arc4_y = teta.map(item => -Math.sin(item))

    var pat1_x: any, pat1_y: any, pat2_x: any, pat2_y: any
    var toReverse: any = []
    if (NP1.rake >= 0) {
        // pat1_x = [plan1_a(:,1); plan1_b(:,1); arc3_x; plan2_a(:,1); plan2_b(:,1); arc1_x];
        pat1_x = []
        plan1_a.elements.forEach((item: any) => pat1_x.push(item[0]))
        plan1_b.elements.forEach((item: any) => pat1_x.push(item[0]))
        arc3_x.forEach(item => pat1_x.push(item))
        plan2_a.elements.forEach((item: any) => pat1_x.push(item[0]))
        plan2_b.elements.forEach((item: any) => pat1_x.push(item[0]))
        arc1_x.forEach(item => pat1_x.push(item))

        // pat1_y = [plan1_a(:,2); plan1_b(:,2); arc3_y; plan2_a(:,2); plan2_b(:,2); arc1_y];
        pat1_y = []
        plan1_a.elements.forEach((item: any) => pat1_y.push(item[1]))
        plan1_b.elements.forEach((item: any) => pat1_y.push(item[1]))
        arc3_y.forEach(item => pat1_y.push(item))
        plan2_a.elements.forEach((item: any) => pat1_y.push(item[1]))
        plan2_b.elements.forEach((item: any) => pat1_y.push(item[1]))
        arc1_y.forEach(item => pat1_y.push(item))

        // pat2_x = [plan2_b(:,1); arc2_x; plan1_b(end:-1:1,1); plan1_a(end:-1:1,1); arc4_x; plan2_a(:,1)];
        pat2_x = []
        plan2_b.elements.forEach((item: any) => pat2_x.push(item[0]))
        arc2_x.forEach(item => pat2_x.push(item))

        plan1_b.elements.forEach((item: any) => toReverse.push(item[0]))
        toReverse.reverse()
        toReverse.forEach((item: any) => pat2_x.push(item))

        toReverse = []
        plan1_a.elements.forEach((item: any) => toReverse.push(item[0]))
        toReverse.reverse()
        toReverse.forEach((item: any) => pat2_x.push(item))

        arc4_x.forEach(item => pat2_x.push(item))

        plan2_a.elements.forEach((item: any) => pat2_x.push(item[0]))

        // pat2_y = [plan2_b(:,2); arc2_y; plan1_b(end:-1:1,2); plan1_a(end:-1:1,2); arc4_y; plan2_a(:,2)];
        pat2_y = []
        plan2_b.elements.forEach((item: any) => pat2_y.push(item[1]))
        arc2_y.forEach(item => pat2_y.push(item))

        toReverse = []
        plan1_b.elements.forEach((item: any) => toReverse.push(item[1]))
        toReverse.reverse()
        toReverse.forEach((item: any) => pat2_y.push(item))

        toReverse = []
        plan1_a.elements.forEach((item: any) => toReverse.push(item[1]))
        toReverse.reverse()
        toReverse.forEach((item: any) => pat2_y.push(item))

        arc4_y.forEach(item => pat2_y.push(item))

        plan2_a.elements.forEach((item: any) => pat2_y.push(item[1]))

        // comp = [pat1_x pat1_y];
        // dilat = [pat2_x pat2_y];
        transposedMatrix = []
        pat1_x.forEach((item: any, idx: number) => transposedMatrix.push([item, pat1_y[idx]]))
        comp = sylvester.Matrix.create(transposedMatrix)

        transposedMatrix = []
        pat2_x.forEach((item: any, idx: number) => transposedMatrix.push([item, pat2_y[idx]]))
        dilat = sylvester.Matrix.create(transposedMatrix)

    } else {
        toReverse = []
        // pat1_x = [plan1_a(:,1); plan1_b(:,1); arc3_x(end:-1:1,1); plan2_b(end:-1:1,1); plan2_a(end:-1:1,1); arc1_x(end:-1:1,1)];
        pat1_x = []
        plan1_a.elements.forEach((item: any) => pat1_x.push(item[0]))
        plan1_b.elements.forEach((item: any) => pat1_x.push(item[0]))

        arc3_x.reverse()
        arc3_x.forEach(item => pat1_x.push(item))

        plan2_b.elements.forEach((item: any) => toReverse.push(item[0]))
        toReverse.reverse()
        toReverse.forEach((item: any) => pat1_x.push(item))

        toReverse = []
        plan2_a.elements.forEach((item: any) => toReverse.push(item[0]))
        toReverse.reverse()
        toReverse.forEach((item: any) => pat1_x.push(item))

        arc1_x.reverse()
        arc1_x.forEach(item => pat1_x.push(item))

        // pat1_y = [plan1_a(:,2); plan1_b(:,2); arc3_y(end:-1:1,1); plan2_b(end:-1:1,2); plan2_a(end:-1:1,2); arc1_y(end:-1:1,1)];
        pat1_y = []

        plan1_a.elements.forEach((item: any) => pat1_y.push(item[1]))
        plan1_b.elements.forEach((item: any) => pat1_y.push(item[1]))

        arc3_y.reverse()
        arc3_y.forEach((item) => pat1_y.push(item))

        toReverse = []
        plan2_b.elements.forEach((item: any) => toReverse.push(item[1]))
        toReverse.reverse()
        toReverse.forEach((item: any) => pat1_y.push(item))

        toReverse = []
        plan2_a.elements.forEach((item: any) => toReverse.push(item[1]))
        toReverse.reverse()
        toReverse.forEach((item: any) => pat1_y.push(item))

        arc1_y.reverse()
        arc1_y.forEach(item => pat1_y.push(item))

        // pat2_x = [plan1_a(:,1); plan1_b(:,1); arc4_x(end:-1:1,1); plan2_a(:,1); plan2_b(:,1); arc2_x(end:-1:1,1)];
        pat2_x = []
        plan1_a.elements.forEach((item: any) => pat2_x.push(item[0]))
        plan1_b.elements.forEach((item: any) => pat2_x.push(item[0]))
        arc4_x.reverse()
        arc4_x.forEach(item => pat2_x.push(item))
        plan2_a.elements.forEach((item: any) => pat2_x.push(item[0]))
        plan2_b.elements.forEach((item: any) => pat2_x.push(item[0]))
        arc2_x.reverse()
        arc2_x.forEach(item => pat2_x.push(item))

        // pat2_y = [plan1_a(:,2); plan1_b(:,2); arc4_y(end:-1:1,1); plan2_a(:,2); plan2_b(:,2); arc2_y(end:-1:1,1)];
        pat2_y = []
        plan1_a.elements.forEach((item: any) => pat2_y.push(item[1]))
        plan1_b.elements.forEach((item: any) => pat2_y.push(item[1]))
        arc4_y.reverse()
        arc4_y.forEach(item => pat2_y.push(item))
        plan2_a.elements.forEach((item: any) => pat2_y.push(item[1]))
        plan2_b.elements.forEach((item: any) => pat2_y.push(item[1]))
        arc2_y.reverse()
        arc2_y.forEach(item => pat2_y.push(item))

        // comp = [pat1_x pat1_y];
        // dilat = [pat2_x pat2_y];
        transposedMatrix = []
        pat1_x.forEach((item: any, idx: number) => transposedMatrix.push([item, pat1_y[idx]]))
        comp = sylvester.Matrix.create(transposedMatrix)

        transposedMatrix = []
        pat2_x.forEach((item: any, idx: number) => transposedMatrix.push([item, pat2_y[idx]]))
        dilat = sylvester.Matrix.create(transposedMatrix);
    }

    return {
        comp: comp,
        dilat: dilat
    }

}

export {doublecouple}
