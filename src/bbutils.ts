import * as _ from 'lodash'

import * as bbmath from './bbmath.ts'
import * as matrices from './matrices.ts'

let numeric = require("numeric");

interface CartesianMomentTensor {
    Mzz: number,
    Mxx: number,
    Myy: number,
    Mxz: number,
    Myz: number,
    Mxy: number
}

interface SphericalMomentTensor {
    Mrr: number,
    Mtt: number,
    Mpp: number,
    Mrt: number,
    Mrp: number,
    Mtp: number
}

/**
 * Convert moment tensor represented by cartesian coordiantes system to
 * moment tensor represented by spherical coordinate system
 *
 * @param mt - moment tensor represented by object 
 * {Mzz: val, Mxx:val, Myy:val, Mxz:val, Myz:val, Mxy:val}
 */
function mtcart2sph(mt: CartesianMomentTensor): SphericalMomentTensor {
    return {
        // cartesian to sphere mapping:
        // cartesian coordinate system [Mzz Mxx Myy Mxz -Myz -Mxy]
        // spherical coordinate system [Mrr Mtt Mpp Mrt  Mrp  Mtp]
        Mrr: mt.Mzz,
        Mtt: mt.Mxx,
        Mpp: mt.Myy,
        Mrt: mt.Mxz,
        Mrp: -mt.Myz,
        Mtp: -mt.Mxy
    }
}

/**
 * Convert moment tensor represented by cartesian coordiantes system to
 * moment tensor represented by spherical coordinate system
 *
 * @param mt - moment tensor represented by object 
 * {Mzz: val, Mxx:val, Myy:val, Mxz:val, Myz:val, Mxy:val}
 */
function mtsph2cart(mt: SphericalMomentTensor): CartesianMomentTensor {
    return {
        // cartesian to sphere mapping:
        // cartesian coordinate system [Mzz Mxx Myy Mxz -Myz -Mxy]
        // spherical coordinate system [Mrr Mtt Mpp Mrt  Mrp  Mtp]
        Mzz: mt.Mrr,
        Mxx: mt.Mtt,
        Myy: mt.Mpp,
        Mxz: mt.Mrt,
        Myz: -mt.Mrp,
        Mxy: -mt.Mtp
    }
}

/**
 * Convert moment tensor represented by spherical coordinate system to the
 * unwrapped symmetric matrix.
 *
 * @param mt - moment tensor represented by object {Mrr: val, Mtt:val, Mpp:val, Mrt:val, Mrp:val, Mtp: val}
 *
 * Format of matrix:
 * [
 *   [Mrr, Mrt, Mrp],
 *   [Mrt, Mtt, Mtp],
 *   [Mrp, Mtp, Mpp]
 * ]
 */
function toMatrix(mt: SphericalMomentTensor) {
    return [
        [mt.Mrr, mt.Mrt, mt.Mrp],
        [mt.Mrt, mt.Mtt, mt.Mtp],
        [mt.Mrp, mt.Mtp, mt.Mpp]
    ]
}

/**
 * Returns object containing eigenvalues and eigenvector for the specified
 * moment tensor
 *
 * @param mt - moment tensor represented by 3x3 matrix in spherical coordinates
 * system
 */
function eigen(mt: number[][]) {
    // check corner case matrices
    var conrnerCaseMatrix = matrices.eig(mt)
    if (conrnerCaseMatrix) {
        return conrnerCaseMatrix
    }

    var eig: { lambda: { x: number[] }, E: { x: number[][] } } = numeric.eig(mt)

    // we need sorted lambda values so keep track of lambda values indexes after sorting
    // and create related eigenvectors matrix
    var sorted = eig.lambda.x.slice(0).sort((a, b) => a - b)
    var indices = sorted.map(item => eig.lambda.x.indexOf(item))

    return {
        eigenvalues: indices.map(i => eig.lambda.x[i]),
        eigenvectors: _.range(indices.length).map(i => indices.map(index => eig.E.x[i][index]))
    }
}

/**
 * Convert cartesian coordinates to spherical coordinates
 *
 * @param x - array of x coordinates
 * @param y - array of y coordinates
 * @param z - arrya of z coordinates
 */
function cart2sph(x: number[], y: number[], z: number[]) {
    var range = _.range(x.length)
    var r = range.map(i => Math.sqrt(x[i] * x[i] + y[i] * y[i] + z[i] * z[i]))
    var elev = range.map(i => Math.atan2(z[i], Math.sqrt(x[i] * x[i] + y[i] * y[i])))
    var az = range.map(i => Math.atan2(y[i], x[i]))
    return {
        r: r,
        elev: elev,
        az: az
    }
}

/**
 * Convert spherical coordinates to cartesian coordinates
 *
 * @param sph - spherical coordinates [azimuth, elevation, radial distance]
 */
function sph2cart(sph: number[]) {
    var [az, elev, r] = sph
    return {
        x: r * Math.cos(elev) * Math.cos(az),
        y: r * Math.cos(elev) * Math.sin(az),
        z: r * Math.sin(elev)
    }
}

/**
 * Convert vector in value-plunge-azimuth to north-east-up
 *
 * @param vpa - [value,plunge,azimuth] array
 */
function vpa2neu(vpa: number[]) {
    var [v, p, a] = vpa
    var enu = sph2cart([Math.PI / 2 - bbmath.D2R * a, -bbmath.D2R * p, v])
    return [enu.y, enu.x, enu.z]
}

/**
 * Get T,P and B axes from moment tensor 
 * 
 * @param mt - moment tensor matrix
 */
function mt2tpb(mt: number[][]) {
    var {
        eigenvalues,
        eigenvectors
    } = eigen(mt)

    var t = [eigenvalues[0], eigenvalues[1], eigenvalues[2]]
    var azelever = cart2sph(eigenvectors[1], numeric.mul(eigenvectors[2], -1), eigenvectors[0])
    var b = azelever.az
    var p = azelever.elev

    // get into degrees
    b = b.map(item => item * bbmath.R2D)

    p = p.map(item => item * bbmath.R2D)

    b = _.zip(p, b).map(pb => pb[0] < 0 ? pb[1] + 180 : pb[1])

    p = p.map(item => Math.abs(item))

    //force azimuth to be in 0-360
    b = b.map(item => item >= 360 ? item - 360 : item)
    b = b.map(item => item < 0 ? item + 360 : item)

    return {
        t: [t[2], p[2], b[2]],
        p: [t[0], p[0], b[0]],
        b: [t[1], p[1], b[1]]
    }
}

/**
 * De-diagonalizes moment tensor into Harvard orientation using eigenvectors in
 * vec
 * @param mt - moment tensor matrix
 * @param vec - eigenvectors
 */
function mt_undiag(mt: number[][], vec: number[][]) {
    mt = numeric.dot(vec, numeric.dot(mt, numeric.transpose(vec)))
    // floating point asymmetry fix
    mt = numeric.add(mt, numeric.transpose(numeric.clone(mt)))
    return mt.map(row => row.map(item => item / 2))
}

/**
 * Decompose moment tensor. Returns a compensated linear
 * vector dipole (clvd) and double couple for the moment tensor
 *
 * @param mt - moment tensor matrix
 */
function mt_decomp(mt: number[][]) {
    var {
        eigenvalues,
        eigenvectors
    } = eigen(mt)

    // first diagonalize (aka principal axis transformation)
    var mt_diag = numeric.diag(eigenvalues)

    var tr = eigenvalues.reduce((acc, cur) => acc + cur) / 3
    var iso = numeric.diag([tr, tr, tr])
    var dev = numeric.sub(mt_diag, iso)

    // only maxdc case
    // identify eigenvalues
    var egn = [dev[0][0], dev[1][1], dev[2][2]]
    var mx = egn.indexOf(Math.max.apply(Math, egn))
    var mn = egn.indexOf(Math.min.apply(Math, egn))

    // best double couple
    var tmp = [0, 0, 0]
    tmp[mx] = (egn[mx] - egn[mn]) / 2
    tmp[mn] = -(egn[mx] - egn[mn]) / 2

    var dblcpl = numeric.diag(tmp)

    // clvd (remainder)
    var clvd = numeric.sub(dev, dblcpl)

    return {
        dblcpl: dblcpl,
        clvd: clvd,
        eigenvectors: eigenvectors
    }
}

/**
 * Converts plane (northing, easting, up) to strike/dip
 *
 * @param neu - [northing, easting, up] array
 */
function norm2strikedip(neu: number[]) {
    var [n, e, u] = neu
    var mod = (x: number, y: number) => x - y * Math.floor(x / y)
    return {
        strike: mod(Math.atan2(-n, e) * bbmath.R2D, 360),
        dip: bbmath.arccosd(u / Math.sqrt(n * n + e * e + u * u))
    }
}

/**
 * Convert normal and slip vectors to strike, dip and rake angles
 * @param normal - normal vector represented as [northing, easting, up]
 * @param slip - slip vector represented as [northing, easting, up]
 */
function normslip2sdr(normal: number[], slip: number[]) {
    // normalize slip vector 
    var smag = Math.sqrt(slip[0] * slip[0] + slip[1] * slip[1] + slip[2] * slip[2])
    slip = slip.map(item => item / smag)

    // get strike & dip from normal    
    var {
        strike,
        dip
    } = norm2strikedip(normal)

    // get rake angle from slip vector and horizontal in-plane vector     
    var comp = bbmath.cosd(strike) * slip[0] + bbmath.sind(strike) * slip[1];
    var rake = bbmath.arccosd(comp);

    return {
        strike: strike,
        dip: dip,
        rake: slip[2] < 0 ? -rake : rake // fix sense of orientation
    }
}

/**
 * Converts principal axes to strike, dip and rake angles
 * @param tension vector - represented as [value, plunge, azimuth]
 * @param compressional vector - represented as [value, plunge, azimuth]
 * @param null vector - represented as [value, plunge, azimuth]
 */
function tpb2sdr(t: number[], p: number[], b: number[]) {
    //force to unit vectors
    t[0] = 1.0
    p[0] = 1.0
    b[0] = 1.0

    t = vpa2neu(t)
    p = vpa2neu(p)
    b = vpa2neu(b)

    var normal: number[] = numeric.add(t, p)
    var slip: number[] = numeric.sub(t, p)

    // precision fix
    // Converting between t and p axes of vertical faults (dip=90deg) and the
    // fault normal can give a slightly downwards normal.  Fix this by
    // adding a small amount to the vertical component.
    if (normal[2] < 0 && normal[2] > -bbmath.EPS) {
        normal[2] = normal[2] + bbmath.EPS
    }

    // if normal is pointing downwards, flip the sign of both the normal and
    // slip vectors to give the appropriate normal and slip vectors (those
    // corresponding to the hanging wall)
    if (normal[2] < 0) {
        normal = normal.map(item => -item)
        slip = slip.map(item => -item)
    }

    // get strike, dip and rake from normal and slip
    var sdr = normslip2sdr(normal, slip)

    // get strike, dip and rake for the swapped normal and slip
    var sdr_swap = normslip2sdr(slip,normal)

    return {
        strike: sdr.strike,
        dip: sdr.dip,
        rake: sdr.rake,
        normalv: normal,
        strike_swap: sdr_swap.strike,
        dip_swap: sdr_swap.dip,
        rake_swap: sdr_swap.rake,
        slipv: slip,
        taxis: [t[1], t[0], t[2]],
        paxis: [p[1], p[0], p[2]],
        baxis: [b[1], b[0], b[2]]
    }
}

/**
 * Converts moment tensor to strike-dip-rake using decomposition
 *
 * @param mt - moment tensor represented by object {Mrr: val, Mtt:val, Mpp:val, 
 * Mrt:val, Mrp:val, Mtp: val}
 */
function mt2sdr(mt: SphericalMomentTensor) {
    var tensormatrix = toMatrix(mt)
    var decomp = mt_decomp(tensormatrix)
    var undiag = mt_undiag(decomp.dblcpl, decomp.eigenvectors)
    var tpb = mt2tpb(undiag)
    var sdr = tpb2sdr(tpb.t, tpb.p, tpb.b)

    return {
        strike: sdr.strike,
        dip: sdr.dip,
        rake: sdr.rake,
        strike_swap: sdr.strike_swap,
        dip_swap: sdr.dip_swap,
        rake_swap: sdr.rake_swap,
        normalv: sdr.normalv,
        slipv: sdr.slipv,
        taxis: sdr.taxis,
        paxis: sdr.paxis,
        baxis: sdr.baxis
    }
}

/**
 * Convert strike, dip and rake to moment tensor represented by matrix in 
 * spherical coordinates system
 *
 * @param sdr - vector [strike, dip, rake]
 */
function sdr2mt(sdr: { strike: number, dip: number, rake: number }): SphericalMomentTensor {
    var { strike, dip, rake} = sdr
    let sq2 = 1/Math.sqrt(2.0)
    return {
        Mrr: sq2 * (bbmath.sind(2 * dip) * bbmath.sind(rake)),
        Mtt: -sq2 * (bbmath.sind(dip) * bbmath.cosd(rake) * bbmath.sind(2 * strike) + bbmath.sind(2 * dip) * bbmath.sind(rake) * Math.pow(bbmath.sind(strike), 2)),
        Mpp: sq2 * (bbmath.sind(dip) * bbmath.cosd(rake) * bbmath.sind(2 * strike) - bbmath.sind(2 * dip) * bbmath.sind(rake) * Math.pow(bbmath.cosd(strike), 2)),
        Mrt: -sq2 * (bbmath.cosd(dip) * bbmath.cosd(rake) * bbmath.cosd(strike) + bbmath.cosd(2 * dip) * bbmath.sind(rake) * bbmath.sind(strike)),
        Mrp: sq2 * (bbmath.cosd(dip) * bbmath.cosd(rake) * bbmath.sind(strike) - bbmath.cosd(2 * dip) * bbmath.sind(rake) * bbmath.cosd(strike)),
        Mtp: -sq2 * (bbmath.sind(dip) * bbmath.cosd(rake) * bbmath.cosd(2 * strike) + 0.5 * bbmath.sind(2 * dip) * bbmath.sind(rake) * bbmath.sind(2 * strike))
    }
}

export { CartesianMomentTensor, SphericalMomentTensor, mt2sdr, mtcart2sph, mtsph2cart, sdr2mt, toMatrix, eigen }
