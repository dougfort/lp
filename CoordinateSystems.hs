{-# OPTIONS -Wall #-}

module CoordinateSystems where

import Diagrams.Backend.Cairo
  ( B,
    renderCairo,
  )
import Diagrams.Prelude
  ( Diagram,
    PolyOrientation (..),
    PolyType (..),
    PolygonOpts (..),
    V2 (..),
    arrowAt,
    black,
    blend,
    deg,
    dims,
    fc,
    lw,
    none,
    p2,
    polygon,
    position,
    r2,
    rad,
    rotate,
    scale,
    sinA,
    white,
    (#),
  )
import MOExamples
  ( Justification (..),
    Table (..),
  )
import Mechanics3D
  ( orient,
    v3FromVec,
  )
import SimpleVec
  ( R,
    Vec,
    iHat,
    jHat,
    kHat,
    magnitude,
    sumV,
    vec,
    xComp,
    yComp,
    zComp,
    zeroV,
    (^/),
  )
import SpatialMath
  ( V3 (..),
  )
import qualified Vis as V

data Position = Cart R R R
  deriving (Show)

type CoordinateSystem = (R, R, R) -> Position

cartesian :: CoordinateSystem
cartesian (x, y, z) = Cart x y z

cart ::
  R -> -- x coordinate
  R -> -- y coordinate
  R -> -- z coordinate
  Position
cart = Cart

cartesianCoordinates :: Position -> (R, R, R)
cartesianCoordinates (Cart x y z) = (x, y, z)

cylindrical :: CoordinateSystem
cylindrical (s, phi, z) = Cart (s * cos phi) (s * sin phi) z

cyl ::
  R -> -- s coordinate
  R -> -- phi coordinate
  R -> -- z coordinate
  Position
cyl s phi z = cylindrical (s, phi, z)

cylindricalCoordinates :: Position -> (R, R, R)
cylindricalCoordinates (Cart x y z) =
  (s, phi, z)
  where
    s = sqrt (x ** 2 + y ** 2)
    phi = atan2 y x

spherical :: CoordinateSystem
spherical (r, theta, phi) =
  Cart (r * sin theta * cos phi) (r * sin theta * sin phi) (r * cos theta)

sph ::
  R -> -- r coordinate
  R -> -- theta coordinate
  R -> -- phi coordinate
  Position
sph r theta phi = spherical (r, theta, phi)

sphericalCoordinates :: Position -> (R, R, R)
sphericalCoordinates (Cart x y z) =
  (r, theta, phi)
  where
    r = sqrt (x ** 2 + y ** 2 + z ** 2)
    theta = atan2 s z
    s = sqrt (x ** 2 + y ** 2)
    phi = atan2 y x

origin :: Position
origin = cart 0 0 0

type Displacement = Vec

displacement ::
  Position -> -- source position
  Position -> -- target position
  Displacement
displacement (Cart x' y' z') (Cart x y z) = vec (x - x') (y - y') (z - z')

shiftPosition :: Displacement -> Position -> Position 
shiftPosition v (Cart x y z) = 
  Cart (x + xComp v) (y + yComp v) (z + zComp v)
