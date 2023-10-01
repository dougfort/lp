{-# OPTIONS -Wall #-}

import CoordinateSystems
  ( VectorField,
    cyl,
    phiHat,
    rVF,
  )
import Geometry
  ( Curve (..),
    Surface (..),
    Volume (..),
  )
import Integrals
  ( crossedLineIntegral,
    curveSample,
    surfaceSample,
    vectorSurfaceIntegral,
    vectorVolumeIntegral,
    volumeSample,
  )
import SimpleVec
  ( R,
    Vec,
    sumV,
    (*^),
    (><),
  )

type Current = R

data CurrentDistribution
  = LineCurrent Current Curve
  | SurfaceCurrent VectorField Surface
  | VolumeCurrent VectorField Volume
  | MultipleCurents [CurrentDistribution]

circularCurrentLoop ::
  R -> -- radius
  R -> -- current
  CurrentDistribution
circularCurrentLoop radius i =
  LineCurrent i (Curve (\phi -> cyl radius phi 0) 0 (2 * pi))

wireSolenoid ::
  R -> -- radius
  R -> -- length
  R -> -- turns/length
  R -> -- current
  CurrentDistribution
wireSolenoid radius len n i =
  LineCurrent
    i
    ( Curve
        (\phi -> cyl radius phi (phi / (2 * pi * n)))
        (-pi * n * len)
        (pi * n * len)
    )

sheetSolenoid ::
  R -> -- radius
  R -> -- length
  R -> -- turns/length
  R -> -- current
  CurrentDistribution
sheetSolenoid radius len n i =
  SurfaceCurrent
    (\r -> (n * i) *^ phiHat r)
    ( Surface
        (uncurry (cyl radius))
        0
        (2 * pi)
        (const $ -len / 2)
        (const $ len / 2)
    )

wireToroid ::
  R -> -- small radius
  R -> -- big radius
  R -> -- number of turns
  R -> -- current
  CurrentDistribution
wireToroid smallR bigR n i =
  let alphi phi = n * phi
      curve phi =
        cyl
          (bigR + smallR * cos (alphi phi))
          phi
          (smallR * sin (alphi phi))
   in LineCurrent i (Curve curve 0 (2 * pi))

magneticDipoleMoment :: CurrentDistribution -> Vec
magneticDipoleMoment (LineCurrent i c) =
  crossedLineIntegral (curveSample 1000) (\r -> 0.5 *^ i *^ rVF r) c
magneticDipoleMoment (SurfaceCurrent k s) =
  vectorSurfaceIntegral (surfaceSample 200) (\r -> 0.5 *^ (rVF r >< k r)) s
magneticDipoleMoment (VolumeCurrent j v) =
  vectorVolumeIntegral (volumeSample 50) (\r -> 0.5 *^ (rVF r >< j r)) v
magneticDipoleMoment (MultipleCurents ds) =
  sumV [magneticDipoleMoment d | d <- ds]
