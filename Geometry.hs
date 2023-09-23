{-# OPTIONS -Wall #-}

module Geometry where

import CoordinateSystems
  ( Position,
    cart,
    cyl,
    cylindrical,
    displacement,
    shiftPosition,
    sph,
    spherical,
  )
import SimpleVec (R, Vec, (*^))

data Curve = Curve
  { curveFunc :: R -> Position,
    startingCurveParam :: R, -- t_a
    endingCurceParam :: R -- t_b
  }

circle2 :: Curve
circle2 = Curve (\t -> cart (2 * cos t) (2 * sin t) 0) 0 (2 * pi)

circle2' :: Curve
circle2' = Curve (\phi -> cyl 2 phi 0) 0 (2 * pi)

unitCircle :: Curve
unitCircle = Curve (\phi -> cyl 1 phi 0) 0 (2 * pi)

straightLine ::
  Position -> -- starting position
  Position -> -- ending position
  Curve -- straight-line curve
straightLine r1 r2 =
  let d = displacement r1 r2
      f t = shiftPosition (t *^ d) r1
   in Curve f 0 1

data Surface = Surface
  { surfaceFunc :: (R, R) -> Position,
    lowerLimit :: R, -- s_l
    upperLimit :: R, -- s_u
    lowerCurve :: R -> R, -- t_l(s)
    upperCurve :: R -> R -- t_u(s)
  }

unitSphere :: Surface
unitSphere =
  Surface
    (\(th, phi) -> sph 1 th phi)
    0
    pi
    (const 0)
    (const $ 2 * pi)
