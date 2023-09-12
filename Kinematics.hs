{-# OPTIONS -Wall #-}

module Kinematics where

import SimpleVec

type Time = R

type Velocity = Vec

type Acceleration = Vec

type VecDerivative = (Time -> Vec) -> Time -> Vec

vecDerivative :: Time -> VecDerivative
vecDerivative dt v t = (v (t + dt / 2) ^-^ v (t - dt / 2)) ^/ dt

velFromPos ::
  R -> -- dt
  (Time -> PosVec) -> -- velocity function
  (Time -> Velocity) -- velocity function
velFromPos = vecDerivative

accFromVel ::
  R -> -- dt
  (Time -> Velocity) -> -- velocity function
  (Time -> Acceleration) -- acceleration function
accFromVel = vecDerivative

positionCV :: PosVec -> Velocity -> Time -> PosVec
positionCV r0 v0 t = v0 ^* t ^+^ r0

velocityCA :: Velocity -> Acceleration -> Time -> Velocity
velocityCA v0 a0 t = a0 ^* t ^+^ v0

positionCA :: PosVec -> Velocity -> Acceleration -> Time -> PosVec
positionCA r0 v0 a0 t = 0.5 *^ t ** 2 *^ a0 ^+^ v0 ^* t ^+^ r0

aParallel :: Vec -> Vec -> Vec
aParallel v a =
  let vHat = v ^/ magnitude v
   in (vHat <.> a) *^ vHat

aPerp :: Vec -> Vec -> Vec
aPerp v a = a ^-^ aParallel v a

speedRateChange :: Vec -> Vec -> R
speedRateChange v a = (v <.> a) / magnitude v

radiusOfCurvature :: Vec -> Vec -> R
radiusOfCurvature v a = (v <.> v) / magnitude (aPerp v a)
