{-# OPTIONS -Wall #-}

import CoordinateSystems
  ( ScalarField,
    VectorField,
    cart,
    magRad,
    rVF,
    shiftPosition,
  )
import qualified Data.Map.Strict as M
import Diagrams.Backend.Cairo (B, renderCairo)
import Diagrams.Prelude
  ( Colour,
    Diagram,
    PolyOrientation (..),
    PolyType (..),
    PolygonOpts (..),
    V2 (..),
    blend,
    deg,
    dims,
    fc,
    lw,
    none,
    p2,
    polygon,
    rad,
    rotate,
    sinA,
    (#),
  )
import qualified Diagrams.Prelude as D
import ElectricField (cSI, mu0)
import SimpleVec
  ( R,
    Vec (..),
    iHat,
    jHat,
    kHat,
    magnitude,
    negateV,
    vec,
    xComp,
    yComp,
    zComp,
    (*^),
    (^+^),
    (^-^),
    (^/),
  )

directionalDerivative :: Vec -> ScalarField -> ScalarField
directionalDerivative d f r =
  (f (shiftPosition (d ^/ 2) r) - f (shiftPosition (negateV d ^/ 2) r))
    / magnitude d

curl :: R -> VectorField -> VectorField
curl a vf r =
  let vx = xComp . vf
      vy = yComp . vf
      vz = zComp . vf
      derivX = directionalDerivative (a *^ iHat)
      derivY = directionalDerivative (a *^ jHat)
      derivZ = directionalDerivative (a *^ kHat)
   in (derivY vz r - derivZ vy r) *^ iHat
        ^+^ (derivZ vx r - derivX vz r) *^ jHat
        ^+^ (derivX vy r - derivY vx r) *^ kHat

type FieldState =
  ( R, -- time t
    VectorField, -- electric Field E
    VectorField -- magnetic Field B
  )

maxwellUpdate ::
  R -> -- dx
  R -> -- dt
  (R -> VectorField) -> -- J
  (FieldState -> FieldState)
maxwellUpdate dx dt j (t, eF, bF) =
  let t' = t + dt
      eF' r = eF r ^+^ cSI ** 2 *^ dt *^ (curl dx bF r ^-^ mu0 *^ j t r)
      bF' r = bF r ^-^ dt *^ curl dx eF r
   in (t', eF', bF')

maxwellEvolve ::
  R -> -- dx
  R -> -- dt
  (R -> VectorField) -> -- J
  (FieldState -> [FieldState])
maxwellEvolve dx dt j = iterate (maxwellUpdate dx dt j)

exLocs, eyLocs, ezLocs, bxLocs, byLocs, bzLocs :: [(Int, Int, Int)]
exLocs = [(nx, ny, nz) | nx <- odds, ny <- evens, nz <- evens]
eyLocs = [(nx, ny, nz) | nx <- evens, ny <- odds, nz <- evens]
ezLocs = [(nx, ny, nz) | nx <- evens, ny <- evens, nz <- odds]
bxLocs = [(nx, ny, nz) | nx <- evens, ny <- odds, nz <- odds]
byLocs = [(nx, ny, nz) | nx <- odds, ny <- evens, nz <- odds]
bzLocs = [(nx, ny, nz) | nx <- odds, ny <- odds, nz <- evens]

spaceStepsCE :: Int
spaceStepsCE = 40

hiEven :: Int
hiEven = 2 * spaceStepsCE

evens :: [Int]
evens = [-hiEven, -hiEven + 2 .. hiEven]

odds :: [Int]
odds = [-hiEven + 1, -hiEven + 3 .. hiEven - 1]
