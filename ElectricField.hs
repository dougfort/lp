{-# OPTIONS -Wall #-}

module ElectricField where

import Charge
  ( Charge,
    ChargeDistribution (..),
    diskCap,
    lineDipole,
    protonOrigin,
    simpleDipole,
  )
import CoordinateSystems
  ( Position,
    ScalarField,
    VectorField,
    addVectorFields,
    cart,
    displacement,
    origin,
    rVF,
    shiftPosition,
    sph,
    vf3D,
    vfGrad,
    vfPNGxy,
  )
import Geometry
  ( Curve (..),
    Surface (..),
    Volume (..),
  )
import SimpleVec
  ( R,
    Vec,
    kHat,
    magnitude,
    sumV,
    vec,
    xComp,
    yComp,
    zComp,
    (*^),
    (<.>),
    (><),
    (^*),
    (^+^),
    (^-^),
    (^/),
  )
