{-# OPTIONS -Wall #-}

module Integrals where

import CoordinateSystems
  ( Position,
    ScalarField,
    displacement,
    origin,
    rVF,
    shiftPosition,
  )
import Geometry
  ( Curve (..),
    Surface (..),
    Volume (..),
  )
import SimpleVec
  ( R,
    Vec,
    magnitude,
    (*^),
    (<.>),
    (><),
    (^+^),
    (^/),
  )

type CurveApprox = Curve -> [(Position, Vec)]

type Segment = (Position, Position)

curveSample :: Int -> Curve -> [(Position, Vec)]
curveSample n c =
  let segCent :: Segment -> Position
      segCent (p1, p2) = shiftPosition ((rVF p1 ^+^ rVF p2) ^/ 2) origin
      segDisp :: Segment -> Vec
      segDisp = uncurry displacement
   in [(segCent seg, segDisp seg) | seg <- segments n c]

segments :: Int -> Curve -> [Segment]
segments n (Curve g a b) =
  let ps = map g $ linSpaced n a b
   in zip ps (tail ps)

linSpaced :: Int -> R -> R -> [R]
linSpaced n x0 x1 =
  take (n + 1) [x0, x0 + dx .. x1]
  where
    dx = (x1 - x0) / fromIntegral n

type ScalerLineIntegral = ScalarField -> Curve -> R

scalarLineIntegral :: CurveApprox -> ScalarField -> Curve -> R
scalarLineIntegral approx f c =
  sum [f r' * magnitude dl' | (r', dl') <- approx c]

type SurfaceApprox = Surface -> [(Position, Vec)]

data Triangle = Tri Position Position Position

triCenter :: Triangle -> Position
triCenter (Tri p1 p2 p3) =
  shiftPosition ((rVF p1 ^+^ rVF p2 ^+^ rVF p3) ^/ 3) origin

triArea :: Triangle -> Vec -- vector area
triArea (Tri p1 p2 p3) = 0.5 *^ (displacement p1 p2 >< displacement p2 p3)

triangles :: Int -> Surface -> [Triangle]
triangles n (Surface g sl su tl tu) =
  let sts = [[(s, t) | t <- linSpaced n (tl s) (tu s)] | s <- linSpaced n sl su]
      stSquares =
        [ ( sts !! j !! k,
            sts !! (j + 1) !! k,
            sts !! (j + 1) !! (k + 1),
            sts !! j !! (k + 1)
          )
          | j <- [0 .. n - 1],
            k <- [0 .. n - 1]
        ]
      twoTriangles (pp1, pp2, pp3, pp4) =
        [Tri (g pp1) (g pp2) (g pp3), Tri (g pp1) (g pp3) (g pp4)]
   in concatMap twoTriangles stSquares

surfaceSample :: Int -> Surface -> [(Position, Vec)]
surfaceSample n s = [(triCenter tri, triArea tri) | tri <- triangles n s]

type ScalarSurfaceIntegral = ScalarField -> Surface -> R

scalarSurfaceIntegral :: SurfaceApprox -> ScalarField -> Surface -> R
scalarSurfaceIntegral approx f s =
  sum [f r' * magnitude da' | (r', da') <- approx s]

type VolumeApprox = Volume -> [(Position, R)]

volumeSample :: Int -> Volume -> [(Position, R)]
volumeSample n v = [(tetCenter tet, tetVolume tet) | tet <- tetrahedrons n v]

data Tet = Tet Position Position Position Position -- tetrahedron

tetCenter :: Tet -> Position
tetCenter (Tet p1 p2 p3 p4) =
  shiftPosition ((rVF p1 ^+^ rVF p2 ^+^ rVF p3 ^+^ rVF p4) ^/ 4) origin

tetVolume :: Tet -> R
tetVolume (Tet p1 p2 p3 p4) =
  abs $ (d1 <.> (d2 >< d3)) / 6
  where
    d1 = displacement p1 p4
    d2 = displacement p2 p4
    d3 = displacement p3 p4

data ParamCube = PC
  { v000 :: (R, R, R),
    v001 :: (R, R, R),
    v010 :: (R, R, R),
    v011 :: (R, R, R),
    v100 :: (R, R, R),
    v101 :: (R, R, R),
    v110 :: (R, R, R),
    v111 :: (R, R, R)
  }

tetrahedrons :: Int -> Volume -> [Tet]
tetrahedrons n (Volume g sl su tl tu ul uu) =
  let stus =
        [ [ [(s, t, u) | u <- linSpaced n (ul s t) (uu s t)]
            | t <- linSpaced n (tl s) (tu s)
          ]
          | s <- linSpaced n sl su
        ]
      stCubes =
        [ PC
            (stus !! j !! k !! l)
            (stus !! j !! k !! (l + 1))
            (stus !! j !! (k + 1) !! l)
            (stus !! j !! (k + 1) !! (l + 1))
            (stus !! (j + 1) !! k !! l)
            (stus !! (j + 1) !! k !! (l + 1))
            (stus !! (j + 1) !! (k + 1) !! l)
            (stus !! (j + 1) !! (k + 1) !! (l + 1))
          | j <- [0 .. n - 1],
            k <- [0 .. n - 1],
            l <- [0 .. n - 1]
        ]
      tets (PC c000 c001 c010 c011 c100 c101 c110 c111) =
        [ Tet (g c000) (g c100) (g c010) (g c001),
          Tet (g c011) (g c111) (g c001) (g c010),
          Tet (g c110) (g c010) (g c100) (g c111),
          Tet (g c101) (g c001) (g c111) (g c100),
          Tet (g c111) (g c100) (g c010) (g c001)
        ]
   in concatMap tets stCubes

type ScalarVolumeIntegral = ScalarField -> Volume -> R

scalarVolumeIntegral :: VolumeApprox -> ScalarField -> Volume -> R
scalarVolumeIntegral approx f vol =
  sum [f r' * dv' | (r', dv') <- approx vol]

vectorLineIntegral = undefined

vectorSurfaceIntegral = undefined

vectorVolumeIntegral = undefined
