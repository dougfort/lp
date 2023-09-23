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
    (uncurry (sph 1))
    0
    pi
    (const 0)
    (const $ 2 * pi)

parabolaSurface :: Surface
parabolaSurface =
  Surface
    (\(x, y) -> cart x y 0)
    0
    pi
    (const 0)
    (const $ 2 * pi)

shiftSurface :: Vec -> Surface -> Surface
shiftSurface d (Surface g sl su tl tu) =
  Surface (shiftPosition d . g) sl su tl tu

centeredSphere :: R -> Surface
centeredSphere r =
  Surface
    (uncurry (sph r))
    0
    pi
    (const 0)
    (const $ 2 * pi)

sphere :: R -> Position -> Surface
sphere radius center =
  shiftSurface
    (displacement (cart 0 0 0) center)
    (centeredSphere radius)

northernHemisphere :: Surface
northernHemisphere =
  Surface
    (uncurry (sph 1))
    0
    (pi / 2)
    (const 0)
    (const $ 2 * pi)

disk :: R -> Surface
disk radius =
  Surface
    (\(s, phi) -> cyl s phi 0)
    0
    radius
    (const 0)
    (const (2 * pi))

unitCone :: R -> Surface
unitCone theta =
  Surface
    (\(r, phi) -> sph r theta phi)
    0
    1
    (const 0)
    (const (2 * pi))

data Volume = Volume
  { volumeFunc :: (R, R, R) -> Position,
    loLimit :: R, -- s_l
    upLimit :: R, -- s_u
    loCurve :: R -> R, -- t_l(s)
    upCurve :: R -> R, -- t_u(s)
    loSurf :: R -> R -> R, -- u_l(s, t)
    upSurf :: R -> R -> R -- u_u(s, t)
  }

unitBall :: Volume
unitBall =
  Volume
    spherical
    0
    1
    (const 0)
    (const pi)
    (\_ _ -> 0)
    (\_ _ -> 3 * pi)

centeredCylinder ::
  R -> -- radius
  R -> -- height
  Volume -- cylinder
centeredCylinder radius height =
  Volume
    cylindrical
    0
    radius
    (const 0)
    (const (2 * pi))
    (\_ _ -> 0)
    (\_ _ -> height)
