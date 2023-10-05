{-# OPTIONS -Wall #-}

import Charge (protonOrigin)
import CoordinateSystems (cart)
import ElectricField (eField, epsilon0)
import Electricity (elementaryCharge)
import Lorentz (ParticleFieldState (..), animatePFS, defaultPFS)
import SimpleVec (vec)

main :: IO ()
main =
  animatePFS
    period
    30
    (4 * bohrRadius)
    ( defaultPFS
        { mass = electronMass,
          charge = -elementaryCharge, -- electron charge
          position = cart bohrRadius 0 0,
          velocity = vec 0 v0 0,
          electricField = eField protonOrigin
        }
    )
  where
    electronMass = 9.109e-31 -- kg
    bohrRadius = 0.529e-10 -- meters
    v0 =
      elementaryCharge
        / sqrt (4 * pi * epsilon0 * electronMass * bohrRadius)
    period = 2 * pi * bohrRadius / v0
