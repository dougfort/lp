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
