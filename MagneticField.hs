{-# OPTIONS -Wall #-}

module MagneticField where

import CoordinateSystems
  ( VectorField,
    addVectorFields,
    cart,
    displacement,
    vfGrad,
  )
import Current
  ( Current,
    CurrentDistribution (..),
    circularCurrentLoop,
    wireSolenoid,
    wireToroid,
  )
import ElectricField
  ( mu0,
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
  ( Vec (..),
    magnitude,
    zComp,
    (*^),
    (><),
    (^/),
  )

bFieldFromLineCurrent ::
  Current -> -- current (in Amps)
  Curve ->
  VectorField -- magnetic field (in Tesla)
bFieldFromLineCurrent i c r =
  let coeff = -mu0 * i / (4 * pi) -- SI units
      integrand r' = d ^/ magnitude d ** 3
        where
          d = displacement r' r
   in coeff *^ crossedLineIntegral (curveSample 1000) integrand c

bFieldFromSurfaceCurrent ::
  VectorField -> -- surface current density (in Amps/m)
  Surface -> -- surface across which the current flows
  VectorField -- magnetic field (in Tesla)
bFieldFromSurfaceCurrent kCurrent s r =
  let coeff = mu0 / (4 * pi) -- SI units
      integrand r' = (kCurrent r' >< d) ^/ magnitude d ** 3
        where
          d = displacement r' r
   in coeff *^ vectorSurfaceIntegral (surfaceSample 200) integrand s

bFieldFromVolumeCurrent ::
  VectorField -> -- volume current density
  Volume -> -- volume throughout which current flows
  VectorField -- magnetic field (in T)
bFieldFromVolumeCurrent j vol r =
  let coeff = mu0 / (4 * pi) -- SI Units
      integrand r' = (j r' >< d) ^/ magnitude d ** 3
        where
          d = displacement r' r
   in coeff *^ vectorVolumeIntegral (volumeSample 50) integrand vol

bField :: CurrentDistribution -> VectorField
bField (LineCurrent i c) = bFieldFromLineCurrent i c
bField (SurfaceCurrent kC s) = bFieldFromSurfaceCurrent kC s
bField (VolumeCurrent j v) = bFieldFromVolumeCurrent j v
bField (MultipleCurents cds) = addVectorFields $ map bField cds

circleB :: VectorField -- magnetic field
circleB = bField $ circularCurrentLoop 0.25 10

bFieldPicSolenoid10 :: IO ()
bFieldPicSolenoid10 =
  vfGrad
    (** 0.2)
    (\(y, z) -> cart 0 (0.02 * y) (0.02 * z))
    (\v -> (yComp v, zComp v))
    "bFieldPicSolenoid10.png"
    20
    (bField $ wireSolenoid 0.01 0.1 100 10)

bFieldPicSolenoid100 :: IO ()
bFieldPicSolenoid100 =
  vfGrad
    (** 0.2)
    (\(y, z) -> cart 0 (0.02 * y) (0.02 * z))
    (\v -> (yComp v, zComp v))
    "bFieldPicSolenoid100.png"
    20
    (bField $ wireSolenoid 0.01 0.1 1000 10)

bFieldWireToroid :: VectorField
bFieldWireToroid =
  bField (wireToroid 0.3 1 50 10)

bFieldPicWireToroid :: IO ()
bFieldPicWireToroid =
  vfGrad
    (** 0.2)
    (\(x, y) -> cart (1.5 * x) (1.5 * y) 0)
    (\v -> (xComp v, yComp v))
    "bFieldPicWireToroid.png"
    20
    bFieldWireToroid
