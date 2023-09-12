{-# OPTIONS -Wall #-}

module Projectile where

import SimpleVec
import Kinematics ( Velocity, PosVec, Time, positionCA )

projectilePos :: PosVec -> Velocity -> Time -> PosVec
projectilePos r0 v0 = positionCA r0 v0 (9.81 *^ negateV kHat)
