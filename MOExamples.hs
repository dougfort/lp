{-# OPTIONS -Wall #-}

module MOExamples where

-- orient

import qualified Graphics.Gloss as G
import Graphics.Gnuplot.Simple
import Mechanics1D
  ( NumericalMethod,
    TimeStep,
    euler,
    rungeKutta4,
  )
import Mechanics3D
  ( HasTime (..),
    ParticleState (..),
    customLabel,
    defaultParticleState,
    disk,
    earthSurfaceGravity,
  )
import MultipleObjects
  ( DMultiParticleState (..),
    Force (..),
    MultiParticleState (..),
    TwoBodyForce,
    billiardForce,
    eulerCromerMPS,
    fixedLinearSpring,
    linearSpring,
    newtonSecondMPS,
    stateMPS,
    updateMPS,
  )
import SimpleVec
  ( R,
    Vec,
    iHat,
    jHat,
    kHat,
    magnitude,
    sumV,
    vec,
    xComp,
    yComp,
    zComp,
    zeroV,
    (*^),
    (^+^),
    (^-^),
  )

twoSpringForces :: [Force]
twoSpringForces =
  [ ExternalForce 0 (fixedLinearSpring 100 0.5 zeroV),
    InternalForce 0 1 (linearSpring 100 0.5),
    ExternalForce 0 earthSurfaceGravity,
    ExternalForce 1 earthSurfaceGravity
  ]

twoSpringInitial :: MultiParticleState
twoSpringInitial =
  MPS
    [ defaultParticleState
        { mass = 2,
          posVec = 0.4 *^ jHat ^-^ 0.3 *^ kHat
        },
      defaultParticleState
        { mass = 3,
          posVec = 0.4 *^ jHat ^-^ 0.8 *^ kHat
        }
    ]

twoSpringUpdate ::
  TimeStep ->
  MultiParticleState -> -- old state
  MultiParticleState -- new state
twoSpringUpdate dt = updateMPS (eulerCromerMPS dt) twoSpringForces
