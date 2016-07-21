/**
 * That's a shame but "numeric" library can't compute eigenvalues and
 * eigenvectors for some "classic" matrices in a reliable way, so we provide
 * precomputed values for these matrices
 */
var cornerCaseMatrices = [{
    matrix: [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0]
    ],
    eigenvalues: [0, 0, 0],
    eigenvectors: [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1]
    ]
}, {
    matrix: [
        [0, 0, 0],
        [0, 0, -1],
        [0, -1, 0]
    ],
    eigenvalues: [-1, 0, 1],
    eigenvectors: [
        [0, 1, 0],
        [-0.7071, 0, -0.7071],
        [-0.7071, 0, 0.7071]
    ]
}, {
    matrix: [
        [0, 0, 0],
        [0, 0, 1],
        [0, 1, 0]
    ],
    eigenvalues: [-1, 0, 1],
    eigenvectors: [
        [0, 1, 0],
        [-0.7071, 0, 0.7071],
        [0.7071, 0, 0.7071]
    ]
}, {
    matrix: [
        [0, -1, 0],
        [-1, 0, 0],
        [0, 0, 0]
    ],
    eigenvalues: [-1, 0, 1],
    eigenvectors: [
        [-0.7071, 0, -0.7071],
        [-0.7071, 0, 0.7071],
        [0, 1, 0]
    ]
}, {
    matrix: [
        [0, 0, 1],
        [0, 0, 0],
        [1, 0, 0]
    ],
    eigenvalues: [-1, 0, 1],
    eigenvectors: [
        [0.7071, 0, 0.7071],
        [0, -1, 0],
        [-0.7071, 0, 0.7071]
    ]
}, {
    matrix: [
        [0, 1, 0],
        [1, 0, 0],
        [0, 0, 0]
    ],
    eigenvalues: [-1, 0, 1],
    eigenvectors: [
        [-0.7071, 0, 0.7071],
        [0.7071, 0, 0.7071],
        [0, 1, 0]
    ]
}, {
    matrix: [
        [0, 0, -1],
        [0, 0, 0],
        [-1, 0, 0]
    ],
    eigenvalues: [-1, 0, 1],
    eigenvectors: [
        [-0.7071, 0, -0.7071],
        [0, 1, 0],
        [-0.7071, 1, 0.7071]
    ]
}];


var isEqual = function(m1:number[][], m2:number[][]) {
    var equal = true;
    var i = 0;
    var j = 0;
    var m1c:string, m2c:string;
    for (; i < 3; i++) {
        j = 0;
        for (; j < 3; j++) {
            m1c = (Math.round(1000 * m1[i][j]) / 1000).toFixed(3);
            m2c = (Math.round(1000 * m2[i][j]) / 1000).toFixed(3);
            if (m1c !== m2c) {
                return false;
            }
        }
    }
    return equal;
};

export function eig(mt:number[][]) {
    for (var i=0 ; i < cornerCaseMatrices.length; i++) {
        let defmat = cornerCaseMatrices[i];
        if (isEqual(defmat.matrix, mt)) {
            return {
                eigenvalues: defmat.eigenvalues.slice(0),
                eigenvectors: defmat.eigenvectors
            };
        }
    }
}
