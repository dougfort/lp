{-# OPTIONS -Wall #-}

import Graphics.Gloss (Picture (..), blue, scale)
import MOExamples
import Mechanics3D (ParticleState (..), simulateGloss)
import MultipleObjects (MultiParticleState (..))
import SimpleVec (iHat, xComp, yComp, zeroV, (*^))

stringPicture :: MultiParticleState -> Picture
stringPicture (MPS sts) =
  let rs = [zeroV] ++ [posVec st | st <- sts] ++ [0.65 *^ iHat]
      xy r = (realToFrac $ xComp r, realToFrac $ yComp r)
      xys = map xy rs
      ppm = 400 -- pixels per meter
   in scale ppm (20 * ppm) $ Color blue $ Line xys

main :: IO ()
main =
  let initialState = stringInitialOvertone 3
   in simulateGloss 0.001 40 initialState stringPicture stringUpdate
