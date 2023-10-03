{-# OPTIONS -Wall #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Lorentz where

import CoordinateSystems
  ( Position (..),
    VectorField,
    addVectorFields,
    cart,
    origin,
    shiftPosition,
    v3FromPos,
    visVec,
  )
import Mechanics1D
  ( Diff (..),
    RealVectorSpace (..),
    rungeKutta4,
  )
import Mechanics3D
  ( HasTime (..),
    simulateVis,
  )
import SimpleVec
  ( R,
    Vec,
    magnitude,
    zeroV,
    (*^),
    (><),
    (^*),
    (^+^),
    (^/),
  )
import qualified Vis as V

data ParticleFieldState = ParticleFieldState
  { mass :: R,
    charge :: R,
    time :: R,
    position :: Position,
    velocity :: Vec,
    electricField :: VectorField,
    magneticField :: VectorField
  }

data DParticleFieldState = DParticleFieldState
  { dmdt :: R,
    dqdt :: R,
    dtdt :: R,
    drdt :: Vec,
    dvdt :: Vec,
    dEdt :: VectorField,
    dBdt :: VectorField
  }

instance RealVectorSpace DParticleFieldState where
  dst1 +++ dst2 =
    DParticleFieldState
      { dmdt = dmdt dst1 + dmdt dst2,
        dqdt = dqdt dst1 + dqdt dst2,
        dtdt = dtdt dst1 + dtdt dst2,
        drdt = drdt dst1 ^+^ drdt dst2,
        dvdt = dvdt dst1 ^+^ dvdt dst2,
        dEdt = addVectorFields [dEdt dst1, dEdt dst2],
        dBdt = addVectorFields [dBdt dst1, dBdt dst2]
      }
  scale w dst =
    DParticleFieldState
      { dmdt = w * dmdt dst,
        dqdt = w * dqdt dst,
        dtdt = w * dtdt dst,
        drdt = w *^ drdt dst,
        dvdt = w *^ dvdt dst,
        dEdt = (w *^) . dEdt dst,
        dBdt = (w *^) . dBdt dst
      }

instance Diff ParticleFieldState DParticleFieldState where
  shift dt dst st =
    ParticleFieldState
      { mass = mass st + dmdt dst * dt,
        charge = charge st + dqdt dst * dt,
        time = time st + dtdt dst * dt,
        position = shiftPosition (drdt dst ^* dt) (position st),
        velocity = velocity st ^+^ dvdt dst ^* dt,
        electricField = \r -> electricField st r ^+^ dEdt dst r ^* dt,
        magneticField = \r -> magneticField st r ^+^ dBdt dst r ^* dt
      }

instance HasTime ParticleFieldState where
  timeOf = time
