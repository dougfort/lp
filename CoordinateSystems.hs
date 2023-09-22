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
    (@@),
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

type ScalarField = Position -> R

addScalarFields :: [ScalarField] -> ScalarField
addScalarFields flds r = sum [fld r | fld <- flds]

type VectorField = Position -> Vec

addVectorFields :: [VectorField] -> VectorField
addVectorFields flds r = sumV [fld r | fld <- flds]

sHat :: VectorField
sHat r =
  vec (cos phi) (sin phi) 0
  where
    (_, phi, _) = cylindricalCoordinates r

phiHat :: VectorField
phiHat r =
  vec (-sin phi) (cos phi) 0
  where
    (_, phi, _) = cylindricalCoordinates r

rHat :: VectorField
rHat rv =
  let d = displacement origin rv
   in if d == zeroV
        then zeroV
        else d ^/ magnitude d

thetaHat :: VectorField
thetaHat r =
  vec
    (cos theta * cos phi)
    (cos theta * sin phi)
    (-sin theta)
  where
    (_, theta, phi) = sphericalCoordinates r

xHat :: VectorField
xHat = const iHat

yHat :: VectorField
yHat = const jHat

zHat :: VectorField
zHat = const kHat

rVF :: VectorField
rVF = displacement origin

sf3d ::
  [Position] -> -- positions to use
  ScalarField -> -- to display
  IO ()
sf3d ps sf =
  V.display whiteBackground $
    orient $
      V.VisObjects
        [ V.Text3d
            (show (round $ sf p :: Int))
            (v3FromPos p)
            V.Fixed9By15
            V.black
          | p <- ps
        ]

v3FromPos :: Position -> V3 R
v3FromPos p =
  V3 x y z
  where
    (x, y, z) = cartesianCoordinates p

whiteBackground :: V.Options
whiteBackground = V.defaultOpts {V.optBackgroundColor = Just V.white}

sfTable ::
  ((R, R) -> Position) ->
  [R] -> -- horizontal
  [R] -> -- vertical
  ScalarField ->
  Table Int
sfTable toPos ss ts sf =
  Table RJ [[round $ sf $ toPos (s, t) | s <- ss] | t <- reverse ts]

vf3D ::
  R -> -- scale factor, vector field units per meter
  [Position] -> -- positions to show the field
  VectorField -> -- vector field to display
  IO ()
vf3D unitsPerMeter ps vf =
  V.display whiteBackground $
    orient $
      V.VisObjects
        [ V.Trans (v3FromPos p) $
            visVec V.black (vf p ^/ unitsPerMeter)
          | p <- ps
        ]

visVec :: V.Color -> Vec -> V.VisObject R
visVec color v =
  let vmag = magnitude v
   in V.Arrow (vmag, 20 * vmag) (v3FromVec v) color

phiHat3D :: IO ()
phiHat3D =
  vf3D
    1
    [ cyl r ph z | r <- [1, 2, 3], ph <- [0, pi / 4 .. 2 * pi], z <- [-2 .. 2]
    ]
    phiHat

vfPNG ::
  ((R, R) -> Position) ->
  (Vec -> (R, R)) ->
  FilePath -> -- file name
  R -> -- scale fascter in units per meter
  [(R, R)] -> -- positions to use
  VectorField -> -- vector field to display
  IO ()
vfPNG toPos fromVec fileName unitsPerMeter pts vf =
  let vf2d = r2 . fromVec . (^/ unitsPerMeter) . vf . toPos
      pic = mconcat [arrowAt (p2 pt) (vf2d pt) | pt <- pts]
   in renderCairo fileName (dims (V2 1024 1024)) pic

vfPNGxy ::
  FilePath -> -- file name
  R -> -- scale fascter in units per meter
  [(R, R)] -> -- positions to use
  VectorField -> -- vector field to display
  IO ()
vfPNGxy = vfPNG (\(x, y) -> cart x y 0) (\v -> (xComp v, yComp v))

vfGrad ::
  (R -> R) ->
  ((R, R) -> Position) ->
  (Vec -> (R, R)) ->
  FilePath ->
  Int -> -- n for n x n
  VectorField ->
  IO ()
vfGrad curve toPos fromVec fileName n vf =
  let step = 2 / fromIntegral n
      xs = [-1 + step / 2, -1 + 3 * step / 2 .. 1 - step / 2]
      pts = [(x, y) | x <- xs, y <- xs]
      array = [(pt, magRad $ fromVec $ vf $ toPos pt) | pt <- pts]
      maxMag = maximum (map (fst . snd) array)
      scaledArrow m th = scale step $ arrowMagRad (curve (m / maxMag)) th
      pic = position [(p2 pt, scaledArrow m th) | (pt, (m, th)) <- array]
   in renderCairo fileName (dims (V2 1024 1024)) pic

magRad :: (R, R) -> (R, R)
magRad (x, y) = (sqrt (x * x + y * y), atan2 y x)

-- magnitude from 0 to 1
-- a picture of an arrow based on a normalized magnitude range
arrowMagRad ::
  R -> -- magnitude
  R -> -- angle in radians, counterclockwise from x axis
  Diagram B
arrowMagRad mag th =
  let r = sinA (15 @@ deg) / sinA (60 @@ deg)
      myType =
        PolyPolar
          [ 120 @@ deg,
            0 @@ deg,
            45 @@ deg,
            30 @@ deg,
            45 @@ deg,
            0 @@ deg,
            120 @@ deg
          ]
          [1, 1, r, 1, 1, r, 1, 1]
      myOpts = PolygonOpts myType NoOrient (p2 (0, 0))
   in scale 0.5 $
        polygon myOpts #lw none # fc (blend mag black white)
          # rotate (th @@ rad)

rVFGrad :: IO ()
rVFGrad =
  vfGrad
    id
    (\(x, y) -> cart x y 0)
    (\v -> (xComp v, yComp v))
    "rVFGrad.png"
    20
    rVF
