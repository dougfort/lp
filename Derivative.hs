{-# OPTIONS -Wall #-}

module Derivative where

import SimpleVec

v1 :: R -> Vec
v1 t = 2 *^ t**2 *^ iHat ^+^ 3 *^ t**3 *^ jHat ^+^ t**4 *^ kHat

xCompFunc :: (R -> Vec) -> R -> R
xCompFunc v t = xComp (v t)

yCompFunc :: (R -> Vec) -> R -> R
yCompFunc v t = yComp (v t)

zCompFunc :: (R -> Vec) -> R -> R
zCompFunc v t = zComp (v t)

type VecDerivative = (R -> Vec) -> R -> Vec

vecDerivative :: R -> VecDerivative
vecDerivative dt v t = (v (t + dt/2) ^-^ v (t - dt/2)) ^/ dt 
