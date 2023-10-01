{-# OPTIONS -Wall #-}

module MagneticField where

import CoordinateSystems
  ( VectorField,
    addVectorFields,
    cart,
    displacement,
    rVF,
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
  ( R,
    Vec (..),
    kHat,
    magnitude,
    zComp,
    (*^),
    (<.>),
    (><),
    (^-^),
    (^/),
  )