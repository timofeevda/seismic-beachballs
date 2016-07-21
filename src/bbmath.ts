export const D2R = Math.PI / 180

export const R2D = 1 / D2R

export const EPS = 2.2204E-016

export const TWO_PI = Math.PI * 2

export const MPI2 = Math.PI / 2

export const MSQRT2 = Math.sqrt(2)

export var sind = (degree:number) => Math.sin(degree * D2R)

export var cosd = (degree:number) => Math.cos(degree * D2R)

export var tand = (degree:number) => Math.tan(degree * D2R)

export var arccosd = (rad: number) => Math.acos(rad) * R2D

export var arctand = (rad: number) => Math.atan(rad) * R2D
