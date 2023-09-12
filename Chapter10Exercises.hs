{-# OPTIONS -Wall #-}

module Chapter10Exercises where

import SimpleVec
import Kinematics


-- Exercise 10.1
v :: R -> Vec
v t = 20 *^ iHat ^-^ (9.8 * t) *^ kHat

r :: R -> Vec
r t = 30 *^ jHat ^+^ (20 *t) *^ iHat ^-^ (4.9 * t**2) *^ kHat

x :: R -> R
x t = iHat <.> r t

-- Exercise 10.2

vecIntegral :: R          -- step size dt 
            -> (R -> Vec) -- vector valued function
            -> R          -- lower limit
            -> R          -- upper limit
            -> Vec        -- result
vecIntegral dt f a b = sumV [f t ^* dt | t <- [a+dt/2, a+3*dt/2 .. b - dt/2]]

-- Exercise 10.3
maxHeight :: PosVec -> Velocity -> R
maxHeight _r0 _v0 = undefined

-- Exercise 10.4
speedCA :: Velocity -> Acceleration -> Time -> R
speedCA v0 a0 t = magnitude $ velocityCA v0 a0 t

-- Exercise 10.5
projectileVel :: Velocity -> Time -> Velocity
projectileVel v0 t = v0 ^+^ (9.81 * t) *^ negateV kHat
