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
  ( CurveApprox,
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
