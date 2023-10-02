{-# OPTIONS -Wall #-}

module Charge where

import CoordinateSystems
  ( Position,
    ScalarField,
    cart,
    displacement,
    origin,
    rVF,
    shiftPosition,
    sph,
  )
import Electricity
  ( elementaryCharge,
  )
import Geometry
  ( Curve (..),
    Surface (..),
    Volume (..),
    disk,
    shiftSurface,
    straightLine,
  )
import Integrals
  ( curveSample,
    scalarLineIntegral,
    scalarSurfaceIntegral,
    scalarVolumeIntegral,
    surfaceSample,
    vectorLineIntegral,
    vectorSurfaceIntegral,
    vectorVolumeIntegral,
    volumeSample,
  )
import SimpleVec
  ( R,
    Vec,
    magnitude,
    negateV,
    sumV,
    vec,
    (*^),
    (<.>),
    (^/),
  )

type Charge = R

data ChargeDistribution
  = PointCharge Charge Position
  | LineCharge ScalarField Curve
  | SurfaceCharge ScalarField Surface
  | VolumeCharge ScalarField Volume
  | MultipleCharges [ChargeDistribution]

protonOrigin :: ChargeDistribution
protonOrigin = PointCharge elementaryCharge origin

chargedLine :: Charge -> R -> ChargeDistribution
chargedLine q len =
  LineCharge (const $ q / len) $
    Curve (cart 0 0) (-len / 2) (len / 2)

chargedBall :: Charge -> R -> ChargeDistribution
chargedBall q radius' =
  VolumeCharge (const $ q / (4 / 3 * pi * radius' ** 3)) $
    Volume
      (\(r, theta, phi) -> sph r theta phi)
      0
      radius'
      (const 0)
      (const pi)
      (\_ _ -> 0)
      (\_ _ -> 2 * pi)

diskCap :: R -> R -> R -> ChargeDistribution
diskCap radius' plateSep sigma =
  MultipleCharges
    [ SurfaceCharge (const sigma) $
        shiftSurface (vec 0 0 (plateSep / 2)) (disk radius'),
      SurfaceCharge (const $ -sigma) $
        shiftSurface (vec 0 0 (-plateSep / 2)) (disk radius')
    ]

totalCharge :: ChargeDistribution -> Charge
totalCharge (PointCharge q _) = q
totalCharge (LineCharge lambda c) =
  scalarLineIntegral (curveSample 1000) lambda c
totalCharge (SurfaceCharge signa s) =
  scalarSurfaceIntegral (surfaceSample 200) signa s
totalCharge (VolumeCharge rho v) =
  scalarVolumeIntegral (volumeSample 50) rho v
totalCharge (MultipleCharges ds) =
  sum [totalCharge d | d <- ds]

simpleDipole ::
  Vec -> -- electric dipole moment
  R -> -- charge separation
  ChargeDistribution
simpleDipole p sep =
  let q = magnitude p / sep
      disp = (sep / 2) *^ (p ^/ magnitude p)
   in MultipleCharges
        [ PointCharge q (shiftPosition disp origin),
          PointCharge (-q) (shiftPosition (negateV disp) origin)
        ]

electricDipoleMoment :: ChargeDistribution -> Vec
electricDipoleMoment (PointCharge q r) =
  q *^ displacement origin r
electricDipoleMoment (LineCharge lambda c) =
  vectorLineIntegral (curveSample 1000) (\r -> lambda r *^ rVF r) c
electricDipoleMoment (SurfaceCharge sigma s) =
  vectorSurfaceIntegral (surfaceSample 200) (\r -> sigma r *^ rVF r) s
electricDipoleMoment (VolumeCharge rho v) =
  vectorVolumeIntegral (volumeSample 50) (\r -> rho r *^ rVF r) v
electricDipoleMoment (MultipleCharges ds) =
  sumV [electricDipoleMoment d | d <- ds]

lineDipole ::
  Vec -> -- dipole moment
  R -> -- charge separation
  ChargeDistribution
lineDipole p sep =
  let disp = (sep / 2) *^ (p ^/ magnitude p)
      curve =
        straightLine
          (shiftPosition (negateV disp) origin)
          (shiftPosition disp origin)
      coeff = 12 / sep ** 3
      lambda r = coeff * (displacement origin r <.> p)
   in LineCharge lambda curve
