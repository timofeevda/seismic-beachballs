import * as bbmath from './bbmath.ts'
import * as bbutils from './bbutils.ts'

let numeric = require("numeric");

interface Ray {
    origin: number[],
    direction: number[]
}

// z axis rotation matrix - 90 degrees clockwise
var ZROT: number[][] = [
    [0, -1, 0],
    [1, 0, 0],
    [0, 0, 1]
]

var identityMatrix = numeric.diag([1, 1, 1])

/**
 * Zero plane normal - used for tessellating polygons intersected by zero plane
 */
var slicingPlaneNormal = [0, 0, 1]

/**
 * Compute cartesian coordinates based on spherical ones
 */
var xyztp = (theta: number, phi: number) => {
    return [
        bbmath.cosd(theta) * bbmath.sind(phi),
        bbmath.sind(theta) * bbmath.sind(phi), bbmath.cosd(phi)
    ]
}

var cross = (a: number[], b: number[]) => {
    return [
        a[1] * b[2] - a[2] * b[1], -(a[0] * b[2] - a[2] * b[0]),
        a[0] * b[1] - a[1] * b[0]
    ]
}

var rot111 = [
    [0, 0, 1],
    [1, 0, 0],
    [0, 1, 0]
]

var rot111dot = [
    [0, 1, 0],
    [0, 0, 1],
    [1, 0, 0]
]

var zref = [
    [1, 0, 0],
    [0, 1, 0],
    [0, 0, -1]
]

var range = (start: number, end: number, step: number) => {
    var range: number[] = []
    while (step > 0 ? end >= start : end <= start) {
        range.push(start)
        start += step
    }
    return range
}


var slice = (start: number, end: number, sliceSize: number) => end <= start ? [] : range(start, end, (end - start) / sliceSize)

var isColored = (vector: number[], lambda: number[]) => numeric.dot(numeric.dot(lambda, vector), vector) > 0 ? true : false

var normalize = (v: number[]) => v.map(i => i / numeric.norm2(v))

// special case for the zero lambda
var th0 = (lambda: number[]) => bbmath.arctand(Math.sqrt(-lambda[0] / lambda[1]))

/**
 * Convert polygon edge to ray
 *
 * @param p1 - origin point
 * @param p2 - destination point
 */
var line2ray = (p1: number[], p2: number[]): Ray => {
    return {
        origin: p1,
        direction: normalize(numeric.sub(p2, p1))
    }
}

/**
 * Find the point where ray is intersected by zero plane
 * @param ray - ray to intersect
 * @param t - distance to plane
 * @returns point of ray and plane intersection
 */
var intersectAt = (ray: Ray, t: number): number[] => numeric.add(numeric.mul(ray.direction, t), ray.origin)

/**
 * Computes distance from plane to the point
 * @param point
 * @returns distance from plane to the point
 */
var distanceToPoint = (point: number[]) => numeric.dot(slicingPlaneNormal, point)

/**
 * Computes distance from ray's origin to the plane
 *
 * @param ray
 * @returns zero if ray is coplanar with plane, null if something strange
 * happened, scalar value representing distance from ray's origin to the
 * plane if everything is OK
 */
function distanceToPlane(ray: Ray) {
    var denominator = numeric.dot(slicingPlaneNormal, ray.direction)
    if (denominator === 0) {
        // ray is coplanar, return 0
        return distanceToPoint(ray.origin) === 0 ? 0 : null
    }

    var t = -(numeric.dot(ray.origin, slicingPlaneNormal)) / denominator
    return t >= 0 ? t : null
}

/**
 * Find point of intersection between ray and zero plane
 *
 * @param ray
 * @returns point of intersection between ray and zero plane
 */
var intersectPlane = (ray: Ray): number[] => intersectAt(ray, distanceToPlane(ray))

function basicPolygons(rangex: number[], rangey: number[]): number[][][] {
    var polygons: number[][][] = []
    for (let i = 0; i < rangex.length - 1; i++) {
        for (let j = 0; j < rangey.length - 1; j++) {
            var vertices: number[][] = []
            vertices.push([rangex[i], rangey[j]])
            vertices.push([rangex[i + 1], rangey[j]])
            vertices.push([rangex[i + 1], rangey[j + 1]])
            vertices.push([rangex[i], rangey[j + 1]])
            polygons.push(vertices)
        }
    }
    return polygons
}

function phiN(theta: number, lambda: number[]) {
    var a = lambda[0] * bbmath.cosd(theta) * bbmath.cosd(theta) + lambda[1] * bbmath.sind(theta) * bbmath.sind(theta)
    var b = a - lambda[2]
    return bbmath.arccosd(Math.sqrt(a / b))
}

function subdivide(poly: number[][], f: (x: number) => number): number[][][] {
    var x1 = poly[0][0]
    var y1 = poly[0][1]
    var x2 = poly[1][0]
    var y2 = poly[1][1]
    var x3 = poly[2][0]
    var y3 = poly[2][1]
    var x4 = poly[3][0]
    var y4 = poly[3][1]

    var f1 = f(x1)
    var f2 = f(x2)
    var f3 = f(x3)
    var f4 = f(x4)

    var slice: number

    if ((y1 <= f1 && f1 <= y4) && (y2 <= f2 && f2 <= y3)) {
        return [
            // lower quadrilateral
            [
                [x1, y1],
                [x2, y2],
                [x2, f2],
                [x1, f1]
            ],
            // upper quadrilateral
            [
                [x4, y4],
                [x3, y3],
                [x2, f2],
                [x1, f1]
            ]
        ]
    } else if ((y1 <= f1 && f1 <= y4) && f2 <= y2) {
        slice = ((y1 - f1) * (x2 - x1)) / (f2 - f1)
        return [
            // upper left triangle
            [
                [x1, f1],
                [x1, y1],
                [x1 + slice, y1],
                [x1, f1]
            ],
            // bottom triangle
            [
                [x4, y4],
                [x2, y2],
                [x3, y3],
                [x4, y4]
            ],
            // remaining 4-sided polygon
            [
                [x4, y4],
                [x2, y2],
                [x1 + slice, y1],
                [x1, f1]
            ]
        ]
    } else if (f1 <= y1 && (y2 <= f2 && f2 <= y3)) {
        slice = ((y1 - f1) * (x2 - x1)) / (f2 - f1)
        return [
            // upper right triangle
            [
                [x1 + slice, y1],
                [x2, y2],
                [x2, f2],
                [x1 + slice, y1]
            ],
            // bottom triangle
            [
                [x4, y4],
                [x1, y1],
                [x3, y3],
                [x4, y4]
            ],
            // remaining 4-sided polygon
            [
                [x1, y1],
                [x3, y3],
                [x2, f2],
                [x1 + slice, y1]
            ]
        ]
    } else if (y4 <= f1 && (y2 <= f2 && f2 <= y3)) {
        slice = ((y3 - f1) * (x2 - x1)) / (f2 - f1)
        return [
            // lower right triangle
            [
                [x1 + slice, y3],
                [x3, y3],
                [x3, f3],
                [x3, f3]
            ],
            // bottom triangle
            [
                [x1, y1],
                [x2, y2],
                [x4, y4],
                [x1, y1]
            ],
            // remaining 4-sided polygon
            [
                [x4, y4],
                [x2, y2],
                [x2, f2],
                [x1 + slice, y3]
            ]
        ]
    } else if ((y1 <= f1 && f1 <= y4) && y3 <= f2) {
        slice = ((y3 - f1) * (x2 - x1)) / (f2 - f1)
        return [
            // lower left triangle
            [
                [x1 + slice, y3],
                [x4, y4],
                [x4, f4],
                [x4, f4]
            ],
            // bottom triangle
            [
                [x1, y1],
                [x2, y2],
                [x3, y3],
                [x1, y1]
            ],
            // remaining 4-sided polygon
            [
                [x1, y1],
                [x3, y3],
                [x1 + slice, y3],
                [x4, f4]
            ]
        ]
    } else if ((f2 <= y2 && y4 <= f4) || (f1 <= y1 && y3 <= f3)) {
        var sliceLeft = ((y1 - f4) * (x2 - x4)) / (f2 - f4)
        var sliceRight = ((y4 - f4) * (x2 - x4)) / (f2 - f4)
        return [
            // left quadrilateral
            [
                [x4 + sliceLeft, y1],
                [x2, y2],
                [x3, y3],
                [x4 + sliceRight, y3]
            ],
            // right quadrilateral
            [
                [x4 + sliceLeft, y1],
                [x1, y1],
                [x4, y4],
                [x4 + sliceRight, y3]
            ]
        ]
    } else {
        // return whole polygon
        return [
            [
                [x1, y1],
                [x2, y2],
                [x3, y3],
                [x4, y4]
            ]
        ]
    }
}

function subdivideForLambda(poly: number[][], lambda: number[]) {
    return subdivide(poly, theta => phiN(theta, lambda))
}

function orangeSlice(theta1: number, theta2: number) {
    var subinterval = Math.min.apply(Math, (range(1, 30, 1)
        .filter(i => (theta2 - theta1) / i <= 10)))

    return basicPolygons(slice(theta1, theta2, subinterval), range(0, 180, 10))
        .map(polygon => polygon.map(vector2 => xyztp(vector2[0], vector2[1])))
}

function subdividePolyList(polyList: number[][][], lambda: number[]) {
    var dividedPolygons: number[][][] = []
    polyList.forEach(polygon => {
        subdivideForLambda(polygon, lambda).forEach(polygon => {
            dividedPolygons.push(polygon)
        })
    })
    return dividedPolygons
}

function polyListForBeachball0(lambda: number[], polyList: number[][][]) {
    var polygons: number[][][] = []
    if ((lambda[0] >= 0 && lambda[1] >= 0 && lambda[2] >= 0) || (lambda[0] <= 0 && lambda[1] <= 0 && lambda[2] <= 0)) {
        polygons = polygons.concat(polyList.map(polygon => {
            return polygon
                .map(point => xyztp(point[0], point[1]))
                .map(point => numeric.dot(zref, point))
        }))

        polygons = polygons.concat(polyList.map(polygon => {
            return polygon
                .map(point => xyztp(point[0], point[1]))
        }))

        return polygons
    } else if (lambda[0] * lambda[1] * lambda[2] === 0) {

        var slice1 = orangeSlice(-th0(lambda), th0(lambda))
        var slice2 = orangeSlice(th0(lambda), -th0(lambda) + 180)
        var slice3 = orangeSlice(-th0(lambda) + 180, th0(lambda) + 180)
        var slice4 = orangeSlice(th0(lambda) + 180, -th0(lambda) + 360)

        polygons = polygons.concat(slice1)
        polygons = polygons.concat(slice2)
        polygons = polygons.concat(slice3)
        polygons = polygons.concat(slice4)

        return polygons
    } else {
        var division = subdividePolyList(polyList, lambda)
        division.forEach(function (polygon) {
            var poly3D = polygon.map(function (point) {
                var point3D = xyztp(point[0], point[1])
                point3D = numeric.dot(zref, point3D)
                return point3D
            })

            polygons.push(poly3D)
        })

        division = subdividePolyList(polyList, lambda)
        division.forEach(function (polygon) {
            var poly3D = polygon.map(function (point) {
                var point3D = xyztp(point[0], point[1])
                return point3D
            })

            polygons.push(poly3D)
        })
        return polygons
    }
}

function cycleMat(vector3: number[]) {
    var [a, b, c] = vector3

    if ((a > 0 && c > 0 && b < 0) || (a < 0 && c < 0 && b > 0) || b === 0) {
        return rot111
    } else if ((b > 0 && c > 0 && a < 0) || (b < 0 && c < 0 && a > 0) || a === 0) {
        return rot111dot
    } else {
        return identityMatrix
    }
}

function polyListForBeachBall(lambda: number[], patternRotation: number[], U: number[][]) {
    // clone lambda
    var lambdaClone = lambda.slice(0)
    var cycleDot = numeric.dot(cycleMat(lambda), lambdaClone)
    var polyList = basicPolygons(range(0, 360, 5), range(0, 90, 5))
    var polys = polyListForBeachball0(cycleDot, polyList)
    var coloredPolys: Array<{ vertices: number[][], compressional: boolean }> = []
    polys.forEach(function (polygon) {
        var mm = numeric.transpose(cycleMat(lambda))
        polygon.forEach(function (point) {
            // rotate with cycle mat
            var p = numeric.dot(mm, point)
            // rotate point
            var rotp = numeric.dot(U, p)
            point[0] = rotp[0]
            point[1] = rotp[1]
            point[2] = rotp[2]
        })

        var xmean = 0
        var ymean = 0
        var zmean = 0
        polygon.forEach(function (point) {
            xmean += point[0]
            ymean += point[1]
            zmean += point[2]
        })
        xmean /= 4
        ymean /= 4
        zmean /= 4

        coloredPolys.push({
            vertices: polygon,
            compressional: isColored([xmean, ymean, zmean], patternRotation)
        })

    })
    return coloredPolys
}

/**
 * Converts strike/dip to fault normal vector
 *
 * @param mt moment tensor represented by object with strike/dip
 * properties
 * @returns fault normal vector in coordinate system where x points east,
 * y points north, z points up
 */
function sd2normal(mt: { dip: number, strike: number }) {
    // compute x,y,z in coordinates system where x points west,
    // y points north and z point up
    var x = -bbmath.sind(mt.dip) * bbmath.sind(mt.strike)
    var y = -bbmath.sind(mt.dip) * bbmath.cosd(mt.strike)
    var z = bbmath.cosd(mt.dip)

    // compute vector position in coordinate system where x points east,
    // y points north, z points up
    return numeric.dot(ZROT, [x, y, z])
}

/**
 * Converts strike/dip/slip to slip vector
 *
 * @param mt - moment tensor represented by object with strike/dip/slip
 * properties
 * @returns slip vector in coordinate system where x points east,
 * y points north, z points up
 */
function sds2slip(mt: { strike: number, dip: number, slip: number }) {
    // compute x,y,z in coordinates system where x points west,
    // y points north and z point up
    var x = bbmath.cosd(mt.slip) * bbmath.cosd(mt.strike) + bbmath.sind(mt.slip) * bbmath.cosd(mt.dip) * bbmath.sind(mt.strike)
    var y = -bbmath.cosd(mt.slip) * bbmath.sind(mt.strike) + bbmath.sind(mt.slip) * bbmath.cosd(mt.dip) * bbmath.cosd(mt.strike)
    var z = bbmath.sind(mt.slip) * bbmath.sind(mt.dip)
    // compute vector position in coordinate system where x points east,
    // y points north, z points up
    return numeric.dot(ZROT, [x, y, z])
}

/**
 * Returns rotation matrix for the beachball polygons. Rotation matrix
 * consists of u1 u2 u3 vectors, which are columns of rotation matrix.
 *
 * Note: u1, u2, u3 are eigenvectors of full moment tensor matrix.
 *
 * @param normal fault normal
 * @param slip slip vector
 * @returns rotation matrix for the beachball polygons
 */
function rotationMatrix(normal: number[], slip: number[]) {
    /**
     * normal = N
     * slip = S
     */

    // u1 = S+N/|S+N|
    var uadd: number[] = numeric.add(slip, normal)
    var u1 = uadd.map(item => item / numeric.norm2(uadd))

    // u3 = S-N/|S-N|
    var usub: number[] = numeric.sub(slip, normal)
    var u3 = usub.map(item => item / numeric.norm2(usub))

    // u2 = u3Xu1 - u2 is a cross product of u3 and u1
    var u2 = cross(u3, u1)

    // transpose to make u1,u2,u3 columns of rotation matrix
    return numeric.transpose([u1, u2, u3])
}

/**
 * Returns rotation matrix for beachball pattern - changes position
 * of colored polygons in beachball
 *
 * @param U rotation matrix consisting of eigenvectors which are columns
 * of rotation matrix. Matrix must be orthonormal, dot(u2,slip) must be
 * equal to 0 (90 degrees between u2 and slip vector).
 *
 * @param lambda vector with moment tensor eigenvalues
 * @returns rotation matrix for beachball pattern - changes position
 * of colored polygons in beachball
 *
 */
function patternRotationMatrix(U: number[][], lambda: number[]) {
    // matrix = dot(dot(U, diag(lambda)),transposed(U)
    return numeric.dot(numeric.dot(U, numeric.diag(lambda)), numeric.transpose(U))
}

/**
 * Converts moment tensor to lambda vector, pattern rotaton matrix and
 * polygons rotation matrix
 *
 * @param mt moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 * @returns object with lambda, patternRotation and rotationMatrix properties
 * filled with corresponding values
 */
function mt2lambda(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    // create symmetric moment tensor matrix based on six components
    var mtMatrix = bbutils.toMatrix(mt)

    // get eigenvalues of moment tensor matrix
    var eigen = bbutils.eigen(mtMatrix)

    // get fault normal and slip vector
    var normal = sd2normal(mt)
    var slip = sds2slip(mt)

    // sort in descending order - required by beachball pattern building
    // algorithm
    var lambda = eigen.eigenvalues.sort((a, b) => b - a)

    // pattern rotation and polygons rotation matrices
    var rotMatrix = rotationMatrix(normal, slip)
    var patRotationMatrix = patternRotationMatrix(rotMatrix, lambda)

    return {
        lambda: lambda,
        patternRotation: patRotationMatrix,
        rotationMatrix: rotMatrix,
        normal: normal,
        slip: slip
    }
}

/**
 * Returns beachball polygons with coloring info
 *
 * @param mt moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns beachball polygons
 */
function beachBall(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    var {lambda, patternRotation, rotationMatrix} = mt2lambda(mt)
    return polyListForBeachBall(lambda, patternRotation, rotationMatrix)
}


/**
 * Split line on first point and intersection point if line is intersected
 * by zero plane
 *
 * @param p1 first point
 * @param p2 second point
 * @param acc vertices array
 */
function splitLine(p1: number[], p2: number[], acc: number[][]) {
    var ray: Ray
    var intersection: number[]
    // check first edge
    if (p1[2] > 0 && p2[2] <= 0) {
        ray = line2ray(p1, p2)
        intersection = intersectPlane(ray)
        acc.push(intersection)
    } else if (p1[2] < 0 && p2[2] >= 0) {
        ray = line2ray(p1, p2)
        intersection = intersectPlane(ray)
        acc.push(p1)
        acc.push(intersection)
    } else {
        // z values are equaly signed and both below zero
        if (p1[2] <= 0 && p2[2] <= 0) {
            acc.push(p1)
        }
    }
}

/**
 * Splits polygon if it is intersected by zero plane
 *
 * @param polygon polygon to split
 * @returns array of polygons - returns array with only one original polygon
 * if it's not intersected by zero plane
 */
function splitPolygon(polygon: { vertices: number[][], compressional: boolean }) {
    var polygons = [polygon]
    var vertices = polygon.vertices
    var p1 = vertices[0]
    var p2 = vertices[1]
    var p3 = vertices[2]
    var p4 = vertices[3]

    var acc: number[][] = []
    splitLine(p1, p2, acc)
    splitLine(p2, p3, acc)
    splitLine(p3, p4, acc)
    splitLine(p4, p1, acc)

    // normalize polygons count
    switch (acc.length) {
        case 3:
            acc.push(acc[0].slice(0))
            polygon.vertices = acc
            break
        case 4:
            polygon.vertices = acc
            break
        case 5:
            // tessellate polygons
            polygon.vertices = acc.slice(0, 4)
            // add second polygon
            polygons.push({
                vertices: [acc[3].slice(0), acc[4].slice(0), acc[0].slice(0), acc[0].slice(0)],
                compressional: polygon.compressional
            })
            break
    }
    return polygons
}

/**
 * Return lower hemisphere of beachball sphere with tessellation of polygons
 * intersected by zero plane
 *
 * @param polygons polygons set representing beachball sphere
 * @returns set of filtered polygons representing lower hemisphere
 */
function filterPolygons(polygons: { vertices: number[][], compressional: boolean }[]) {
    var filteredPolygons: { vertices: number[][], compressional: boolean }[] = []

    polygons.forEach(polygon => {
        // check if all vertices are below zero
        var lower = polygon.vertices[0][2] <= 0 && polygon.vertices[1][2] <= 0 && polygon.vertices[2][2] <= 0 && polygon.vertices[3][2] <= 0

        // check if all vertices are above zero
        var higher = polygon.vertices[0][2] >= 0 && polygon.vertices[1][2] >= 0 && polygon.vertices[2][2] >= 0 && polygon.vertices[3][2] >= 0

        // add polygon if all vertices are below zero
        if (lower) {
            filteredPolygons.push(polygon)
        } else if (!higher) {
            // this means that polygon is intersected by zero plane, so we
            // need to split polygons intersected by zero plane (kind
            // of tessellation to improve rendering quality)
            splitPolygon(polygon).forEach(polygon => {
                filteredPolygons.push(polygon)
            })
        }
    })

    return filteredPolygons
}

/**
 * Return polygons representing lower hemisphere of beachball
 *
 * @param mt moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns polygons set representing lower hemisphere of beachball
 */
function lowerHemisphere(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    return filterPolygons(beachBall(mt))
}

/**
 * Return polygons representing upper hemisphere of beachball
 *
 * @param mt moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns upper hemisphere of beachball
 */
function upperHemisphere(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    var polygons = beachBall(mt)

    // negate vertices to treat upper hemisphere in the same way
    // as lower one
    polygons.forEach(polygon => {
        polygon.vertices.forEach(point => {
            point[2] = -point[2]
        })
    })

    var filteredPolygons = filterPolygons(polygons)

    // negate back to transform lower hemisphere to upper hemisphere
    filteredPolygons.forEach(polygon => {
        polygon.vertices.forEach(point => {
            point[2] = -point[2]
        })
    })

    return filteredPolygons
}

/**
 * Modifies vertex applying Wulff Net projection
 *
 * @param vertex vertex to modify
 */
function wulffNetProjection(vertex: number[]) {
    vertex[0] = vertex[0] / (1 - vertex[2])
    vertex[1] = vertex[1] / (1 - vertex[2])
    vertex[2] = 0
}

/**
 * Modifies vertex applying Equal Area Net projection
 *
 * @param vertex vertex to modify
 */
function equalAreaNetProjection(vertex: number[]) {
    vertex[0] = vertex[0] * Math.sqrt(1 / (1 - vertex[2]))
    vertex[1] = vertex[1] * Math.sqrt(1 / (1 - vertex[2]))
    vertex[2] = 0
}

/**
 * Modifies vertex applying orthographic projection
 *
 * @param vertex vertex to modify
 */
function orthographicProjection(vertex: number[]) {
    vertex[2] = 0
}

/**
 * Return polygons representing lower hemisphere of beachball. Polygons
 * are project on horizontal plane - z values are zeroed
 *
 * @param mt moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns lower hemisphere of beachball with polygons projected on
 * horizontal plane with orthographic projection
 */
function lowerHemisphereOrthographic(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    var lHemisphere = lowerHemisphere(mt)
    lHemisphere.forEach(polygon => {
        polygon.vertices.forEach(vertex => {
            orthographicProjection(vertex)
        })
    })
    return lHemisphere
}

/**
 * Return polygons representing lower hemisphere of beachball. Polygons
 * are project on horizontal plane using Wulff Net stereographic projection
 *
 * @param mt moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns lower hemisphere of beachball with polygons projected on
 * horizontal plane with Wulff Net stereographic projection
 */
function lowerHemisphereWulffNet(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    var lHemisphere = lowerHemisphere(mt)
    lHemisphere.forEach(polygon => {
        polygon.vertices.forEach(vertex => {
            wulffNetProjection(vertex)
        })
    })
    return lHemisphere
}

/**
 * Return polygons representing lower hemisphere of beachball. Polygons
 * are project on horizontal plane using Equal Area Net stereographic projection
 *
 * @param mt - moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns lower hemisphere of beachball with polygons projected on
 * horizontal plane with Equal Area Net stereographic projection
 */
function lowerHemisphereEqualAreaNet(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    var lHemisphere = lowerHemisphere(mt)
    lHemisphere.forEach(polygon => {
        polygon.vertices.forEach(vertex => {
            equalAreaNetProjection(vertex)
        })
    })
    return lHemisphere
}

/**
 * Return polygons representing upper hemisphere of beachball. Polygons
 * are project on horizontal plane - z values are zeroed
 *
 * @param mt - moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns upper hemisphere of beachball with polygons projected on
 * horizontal plane with orthographic projection
 */
function upperHemisphereOrthographic(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    var uHemisphere = upperHemisphere(mt)
    // zero all z values
    uHemisphere.forEach(polygon => {
        polygon.vertices.forEach(vertex => {
            vertex[2] = -vertex[2]
            orthographicProjection(vertex)
        })
    })
    return uHemisphere
}

/**
 * Return polygons representing upper hemisphere of beachball. Polygons
 * are project on horizontal plane using Wulff Net stereographic projection
 *
 * @param mt moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns upper hemisphere of beachball with polygons projected on
 * horizontal plane with Wulff Net stereographic projection
 */
function upperHemisphereWulffNet(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    var uHemisphere = upperHemisphere(mt)
    // zero all z values
    uHemisphere.forEach(polygon => {
        polygon.vertices.forEach(vertex => {
            vertex[2] = -vertex[2]
            wulffNetProjection(vertex)
        })
    })
    return uHemisphere
}

/**
 * Return polygons representing lower hemisphere of beachball. Polygons
 * are project on horizontal plane using Equal Area Net stereographic projection
 *
 * @param mt moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns lower hemisphere of beachball with polygons projected on
 * horizontal plane with Equal Area Net stereographic projection
 */
function upperHemisphereEqualAreaNet(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    var uHemisphere = upperHemisphere(mt)
    // zero all z values
    uHemisphere.forEach(polygon => {
        polygon.vertices.forEach(vertex => {
            vertex[2] = -vertex[2]
            equalAreaNetProjection(vertex)
        })
    })
    return uHemisphere
}

/**
 * Return fault normal and slip vectors
 *
 * @param mt moment tensor object with 6 moment tensor components and
 * strike/dip/slip as properties
 *
 * @returns object with normal and slip vectors
 */
function normalslip(mt: bbutils.SphericalMomentTensor & { strike: number, dip: number, slip: number }) {
    var {normal, slip} = mt2lambda(mt)
    return { normal, slip }
}

export {
beachBall,
lowerHemisphere,
upperHemisphere,
lowerHemisphereEqualAreaNet,
lowerHemisphereWulffNet,
lowerHemisphereOrthographic,
upperHemisphereEqualAreaNet,
upperHemisphereWulffNet,
upperHemisphereOrthographic,
normalslip
}
