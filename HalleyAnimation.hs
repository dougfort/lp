{-# OPTIONS -WALL #-}

import Graphics.Gloss (Picture (..), pictures, red, translate, yellow)
import Mechanics3D
  ( ParticleState (..),
    disk,
    halleyInitial,
    halleyUpdate,
    simulateGloss,
  )
import SimpleVec (xComp, yComp)

diskComet :: Picture
diskComet = Color red (disk 10)

diskSun :: Picture
diskSun = Color yellow (disk 20)

halleyPicture :: ParticleState -> Picture
halleyPicture (ParticleState _m _q _t r _v) =
  pictures [diskSun, translate xPixels yPixels diskComet] 
    where 
      pixelsPerMeter = 1e-10
      xPixels = pixelsPerMeter * realToFrac (xComp r)
      yPixels = pixelsPerMeter * realToFrac (yComp r)

main :: IO ()
main = simulateGloss (365.25 * 24 * 60 * 60) 400 
  halleyInitial halleyPicture halleyUpdate
