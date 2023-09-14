{-# OPTIONS -Wall #-}

import MOExamples (twoSpringInitial, twoSpringUpdate)
import Mechanics3D (posVec, simulateVis, v3FromVec)
import MultipleObjects (MultiParticleState (..))
import SimpleVec (R, zeroV)
import Vis (Flavour (..), VisObject (..), blue, green, red)

main :: IO ()
main = simulateVis 1 20 twoSpringInitial twoSpringVisObject twoSpringUpdate

twoSpringVisObject :: MultiParticleState -> VisObject R
twoSpringVisObject (MPS sts) =
  let r0 = posVec (head sts)
      r1 = posVec (sts !! 1)
      springsObj =
        Line
          Nothing
          [ v3FromVec zeroV,
            v3FromVec r0,
            v3FromVec r1
          ]
          blue
      objs =
        [ Trans (v3FromVec r0) (Sphere 0.1 Solid red),
          Trans (v3FromVec r1) (Sphere 0.1 Solid green),
          springsObj
        ]
      vpm = 1 -- Vis units per meter
   in Scale (vpm, vpm, vpm) $ VisObjects objs
