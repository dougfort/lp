{-# OPTIONS -Wall #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# OPTIONS_GHC -Wno-unrecognised-pragmas #-}

{-# HLINT ignore "Use newtype instead of data" #-}
{-# HLINT ignore "Use zipWith" #-}

module MultipleObjects where

import Data.List (zipWith)
import GHC.Base (undefined)
import Mechanics1D
  ( Diff (..),
    Mass,
    NumericalMethod,
    RealVectorSpace (..),
    TimeStep,
    euler,
  )
import Mechanics3D
  ( DParticleState (..),
    HasTime (..),
    OneBodyForce,
    ParticleState (..),
    defaultParticleState,
    newtonSecondPS,
  )
import SimpleVec
  ( R,
    Vec,
    magnitude,
    zeroV,
    (*^),
    (<.>),
    (><),
    (^*),
    (^+^),
    (^-^),
    (^/),
  )

type TwoBodyForce =
  ParticleState -> -- for is produced BY particle with this state
  ParticleState -> -- force acts ON particle with this state
  ForceVector

type ForceVector = Vec

oneFromTwo ::
  ParticleState -> -- state of the particle PRODUCING the force
  TwoBodyForce ->
  OneBodyForce
oneFromTwo stBy f = f stBy

gravityMagnitude :: Mass -> Mass -> R -> R
gravityMagnitude m1 m2 r =
  let gg = 6.67408e-11 -- N m^2 / kg^2
   in gg * m1 * m2 / r ** 2

universalGravity :: TwoBodyForce
universalGravity st1 st2 =
  let gg = 6.67408e-11 -- N m^2 / kg^2
      m1 = mass st1
      m2 = mass st2
      r1 = posVec st1
      r2 = posVec st2
      r21 = r2 ^-^ r1
   in (-gg) *^ m1 *^ m2 *^ r21 ^/ magnitude r21 ** 3

constantRepulsiveForce :: R -> TwoBodyForce
constantRepulsiveForce force st1 st2 =
  let r1 = posVec st1
      r2 = posVec st2
      r21 = r2 ^-^ r1
   in force *^ r21 ^/ magnitude r21

linearSpring ::
  R -> -- spring constant
  R -> -- equilibrium length
  TwoBodyForce
linearSpring k re st1 st2 =
  let r1 = posVec st1
      r2 = posVec st2
      r21 = r2 ^-^ r1
      r21mag = magnitude r21
   in (-k) *^ (r21mag - re) *^ r21 ^/ r21mag

fixedLinearSpring :: R -> R -> Vec -> OneBodyForce
fixedLinearSpring k re r1 =
  oneFromTwo (defaultParticleState {posVec = r1}) (linearSpring k re)

centralForce :: (R -> R) -> TwoBodyForce
centralForce f st1 st2 =
  let r1 = posVec st1
      r2 = posVec st2
      r21 = r2 ^-^ r1
      r21mag = magnitude r21
   in f r21mag *^ r21 ^/ r21mag

linearSpringCentral ::
  R -> -- spring constant
  R -> -- equilibrium length
  TwoBodyForce
linearSpringCentral k re = centralForce (\r -> -k * (r - re))

universalGravityCentral :: TwoBodyForce
universalGravityCentral st1 st2 =
  let gg = 6.67408e-11 -- N m^2 / kg^2
      m1 = mass st1
      m2 = mass st2
   in centralForce (\r -> (-gg) * m1 * m2 / r ** 2) st1 st2

billiardForce ::
  R -> -- spring constant
  R -> -- threshold center separation
  TwoBodyForce
billiardForce k re =
  centralForce $ \r ->
    if r >= re
      then 0
      else (-k * (r - re))

data Force
  = ExternalForce Int OneBodyForce
  | InternalForce Int Int TwoBodyForce

data MultiParticleState = MPS {particleStates :: [ParticleState]} deriving (Show)

instance HasTime MultiParticleState where
  timeOf (MPS sts) = time (head sts)

data DMultiParticleState = DMPS [DParticleState] deriving (Show)

instance RealVectorSpace DMultiParticleState where
  DMPS dsts1 +++ DMPS dsts2 = DMPS $ zipWith (+++) dsts1 dsts2
  scale w (DMPS dsts) = DMPS $ map (scale w) dsts

instance Diff MultiParticleState DMultiParticleState where
  shift dt (DMPS dsts) (MPS sts) = MPS $ zipWith (shift dt) dsts sts

newtonSecondMPS ::
  [Force] ->
  (MultiParticleState -> DMultiParticleState) -- differential equation
newtonSecondMPS fs mpst@(MPS sts) =
  let deriv (n, st) = newtonSecondPS (forcesOn n mpst fs) st
   in DMPS $ map deriv (zip [0 ..] sts)

forcesOn :: Int -> MultiParticleState -> [Force] -> [OneBodyForce]
forcesOn n mpst = map (forceOn n mpst)

forceOn :: Int -> MultiParticleState -> Force -> OneBodyForce
forceOn n _ (ExternalForce n0 fOneBody)
  | n == n0 = fOneBody
  | otherwise = const zeroV
forceOn n (MPS sts) (InternalForce n0 n1 fTwoBody)
  | n == n0 = oneFromTwo (sts !! n1) fTwoBody -- n1 acts on n0
  | n == n1 = oneFromTwo (sts !! n0) fTwoBody -- n0 acts on n1
  | otherwise = const zeroV

eulerCromerMPS ::
  TimeStep -> -- dt for stepping
  NumericalMethod MultiParticleState DMultiParticleState
eulerCromerMPS dt deriv mpst0 =
  let mpst1 = euler dt deriv mpst0
      sts0 = particleStates mpst0
      sts1 = particleStates mpst1
   in -- now update positions
      MPS $
        [ st1 {posVec = posVec st0 ^+^ velocity st1 ^* dt}
          | (st0, st1) <- zip sts0 sts1
        ]

updateMPS ::
  NumericalMethod MultiParticleState DMultiParticleState ->
  [Force] ->
  (MultiParticleState -> MultiParticleState)
updateMPS method = method . newtonSecondMPS

stateMPS ::
  NumericalMethod MultiParticleState DMultiParticleState ->
  [Force] ->
  (MultiParticleState -> [MultiParticleState])
stateMPS method = iterate . method . newtonSecondMPS
