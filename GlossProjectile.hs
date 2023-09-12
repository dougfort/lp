{-# OPTIONS -Wall #-}

import Graphics.Gloss (Picture (..), red, scale, translate)
import Mechanics3D
  ( ParticleState (..),
    disk,
    projectileInitial,
    projectileUpdate,
    simulateGloss,
  )
import SimpleVec (yComp, zComp)
import System.Environment (getArgs)

projectilePicture :: ParticleState -> Picture
projectilePicture (ParticleState _m _q _t r _v) =
  scale 0.2 0.2 $ translate yFloat zFloat redDisk
  where
    yFloat = realToFrac (yComp r)
    zFloat = realToFrac (zComp r)
    redDisk :: Picture
    redDisk = Color red (disk 50)

mainWithArgs :: [String] -> IO ()
mainWithArgs args =
  simulateGloss 3 20 (projectileInitial args) projectilePicture projectileUpdate

main :: IO ()
main = getArgs >>= mainWithArgs
