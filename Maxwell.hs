{-# OPTIONS -Wall #-}
{-# OPTIONS_GHC -Wno-unrecognised-pragmas #-}

{-# HLINT ignore "Replace case with fromMaybe" #-}

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

-- Finite Difference Time Domain
data StateFDTD = StateFDTD
  { timeFDTD :: R,
    stepX :: R,
    stepY :: R,
    stepZ :: R,
    eField :: M.Map (Int, Int, Int) R,
    bField :: M.Map (Int, Int, Int) R
  }
  deriving (Show)

initialStateFDTD :: R -> StateFDTD
initialStateFDTD spatialStep =
  StateFDTD
    { timeFDTD = 0,
      stepX = spatialStep,
      stepY = spatialStep,
      stepZ = spatialStep,
      eField = M.fromList [(loc, 0) | loc <- exLocs ++ eyLocs ++ ezLocs],
      bField = M.fromList [(loc, 0) | loc <- bxLocs ++ byLocs ++ bzLocs]
    }

lookupAZ :: Ord k => k -> M.Map k R -> R
lookupAZ key m =
  case M.lookup key m of
    Just x -> x
    Nothing -> 0

partialX, partialY, partialZ :: R -> M.Map (Int, Int, Int) R -> (Int, Int, Int) -> R
partialX dx m (i, j, k) = (lookupAZ (i + 1, j, k) m - lookupAZ (i - 1, j, k) m) / dx
partialY dx m (i, j, k) = (lookupAZ (i, j + 1, k) m - lookupAZ (i, j - 1, k) m) / dx
partialZ dx m (i, j, k) = (lookupAZ (i, j, k + 1) m - lookupAZ (i, j, k - 1) m) / dx

curlEx, curlEy, curlEz, curlBx, curlBy, curlBz :: StateFDTD -> (Int, Int, Int) -> R
curlEx (StateFDTD _ _ dy dz e _) loc = partialY dy e loc - partialZ dz e loc
curlEy (StateFDTD _ dx _ dz e _) loc = partialZ dz e loc - partialX dx e loc
curlEz (StateFDTD _ dx dy _ e _) loc = partialX dx e loc - partialY dy e loc
curlBx (StateFDTD _ _ dy dz _ b) loc = partialY dy b loc - partialZ dz b loc
curlBy (StateFDTD _ dx _ dz _ b) loc = partialZ dz b loc - partialX dx b loc
curlBz (StateFDTD _ dx dy _ _ b) loc = partialX dx b loc - partialY dy b loc

stateUpdate :: 
  R -> -- dt
  (R -> VectorField) ->  -- current density
  (StateFDTD -> StateFDTD)
stateUpdate dt j st0@(StateFDTD t _dx _dy _dz _e _b) =
  let st1 = updateE dy (j t) st0 
      st2 = updateB dt st1 
   in st2 
    