{-# OPTIONS -Wall #-}

module Electricity where

import Graphics.Gnuplot.Simple
  ( Attribute (..),
    plotPaths,
  )
import Mechanics3D
  ( ParticleState (..),
    defaultParticleState,
  )
import MultipleObjects
  ( Force (..),
    MultiParticleState (..),
    TwoBodyForce,
    centralForce,
    eulerCromerMPS,
    statesMPS,
  )
import SimpleVec
  ( R,
    Vec (..),
    iHat,
    (*^),
  )

type Charge = R

elementaryCharge :: Charge
elementaryCharge = 1.602176634e-19 -- in Coulombs

coulombMagnitude :: Charge -> Charge -> R -> R
coulombMagnitude q1 q2 r =
  let k = 9e9 -- N m^2 / C^2
   in k * abs (q1 * q2) / r ** 2

coulombForce :: TwoBodyForce
coulombForce st1 st2 =
  let k = 9e9 -- N m^2 / C^2
      q1 = charge st1
      q2 = charge st2
   in centralForce (\r -> k * q1 * q2 / r ** 2) st1 st2

twoProtonStates ::
  R -> -- time step
  MultiParticleState -> -- initial 2 particle state
  [MultiParticleState] -- infinite list of states
twoProtonStates dt =
  statesMPS (eulerCromerMPS dt) [InternalForce 1 0 coulombForce]

initialTwoProtonState ::
  R -> -- initial separation
  MultiParticleState
initialTwoProtonState d =
  let protonMass = 1.673e-27 -- in kg
   in MPS
        [ defaultParticleState
            { mass = protonMass,
              charge = elementaryCharge,
              posVec = (-d / 2) *^ iHat
            },
          defaultParticleState
            { mass = protonMass,
              charge = elementaryCharge,
              posVec = (d / 2) *^ iHat
            }
        ]

oneProtonVelocity ::
  R -> -- dt
  R -> -- starting separation
  [(R, R)] -- (time, velocity) pairs
oneProtonVelocity dt d =
  let state0 = initialTwoProtonState d
   in [ (time st2, xComp $ velocity st2)
        | MPS [_, st2] <- twoProtonStates dt state0
      ]

tvPairs :: [(R, R)]
tvPairs = takeWhile (\(t, _) -> t <= 2e-2) $ oneProtonVelocity 1e-5 1e-2

velocityPlot :: IO ()
velocityPlot =
  plotPaths
    [ Title "Two protons released from 1 cm",
      XLabel "Time (s)",
      YLabel "Proton velocity (m/s)",
      PNG "protons.png",
      Key Nothing
    ]
    [ tvPairs,
      [(t, 1379 * t) | t <- [0, 1e-5 .. 4e-3]],
      [(t, 3.71) | t <- [0, 1e-3 .. 2e-2]]
    ]
