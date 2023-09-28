{-# OPTIONS -Wall #-}

module ElectricField where

import Charge
  ( Charge,
    ChargeDistribution (..),
    diskCap,
    lineDipole,
    protonOrigin,
    simpleDipole,
  )
import CoordinateSystems
  ( Position,
    ScalarField,
    VectorField,
    addVectorFields,
    cart,
    displacement,
    origin,
    rVF,
    shiftPosition,
    sph,
    vf3D,
    vfGrad,
    vfPNGxy,
  )
import Diagrams (r2)
import Geometry
  ( Curve (..),
    Surface (..),
    Volume (..),
  )
import Integrals
  ( curveSample,
    surfaceSample,
    vectorLineIntegral,
    vectorSurfaceIntegral,
  )
import SimpleVec
  ( R,
    Vec,
    kHat,
    magnitude,
    sumV,
    vec,
    xComp,
    yComp,
    zComp,
    (*^),
    (<.>),
    (><),
    (^*),
    (^+^),
    (^-^),
    (^/),
  )

epsilon0 :: R
epsilon0 = 1 / (mu0 * cSI ** 2)

cSI :: R
cSI = 299792458 -- m/s

mu0 :: R
mu0 = 4e-7 * pi -- N/A*2

eFieldFromPointCharge ::
  Charge -> -- in Coulombs
  Position -> -- of point charge (in m)
  VectorField -- electric field (in V/m)
eFieldFromPointCharge q1 r1 r =
  let k = 1 / (4 * pi * epsilon0)
      d = displacement r1 r
   in (k * q1) *^ d ^/ magnitude d ** 3

eFieldFromLineCharge ::
  ScalarField -> -- linear charge density lambda
  Curve -> -- geometry of the line charge
  VectorField -- electric field (in V/m)
eFieldFromLineCharge lambda c r =
  let k = 1 / (4 * pi * epsilon0)
      integrand r' = lambda r' *^ d ^/ magnitude d ** 3
        where
          d = displacement r' r
   in k *^ vectorLineIntegral (curveSample 1000) integrand c

eFieldFromSurfaceCharge ::
  ScalarField -> -- surface charge density sigma
  Surface -> -- geometry of the surface charge
  VectorField -- electric field (in V/m)
eFieldFromSurfaceCharge sigma s r =
  let k = 1 / (4 * pi * epsilon0)
      integrand r' = sigma r' *^ d ^/ magnitude d ** 3
        where
          d = displacement r' r
   in k *^ vectorSurfaceIntegral (surfaceSample 200) integrand s

eField :: ChargeDistribution -> VectorField
eField (PointCharge q r) = eFieldFromPointCharge q r
eField (LineCharge lam c) = eFieldFromLineCharge lam c
eField (SurfaceCharge sig s) = eFieldFromSurfaceCharge sig s
eField (MultipleCharges cds) = addVectorFields $ map eField cds -- superposition

eFieldPicProton2D :: IO ()
eFieldPicProton2D =
  vfPNGxy "eFieldPicProton2D.png" 3e-9 pts (eField protonOrigin)
  where
    pts =
      [ (r * cos th, r * sin th) | r <- [1, 1.5, 2], th <- [0, pi / 4 .. 2 * pi]
      ]

eFieldPicProtonGrad :: IO ()
eFieldPicProtonGrad =
  vfGrad
    (** 0.2)
    (\(x, y) -> cart x y 0)
    (\v -> (xComp v, yComp v))
    "eFieldPicProtonGrad.png"
    20
    (eField protonOrigin)

eFieldPicProton3D :: IO ()
eFieldPicProton3D =
  vf3D
    4e-9
    [ sph r th ph | r <- [1, 1.5, 2], th <- [0, pi / 4 .. pi], ph <- [0, pi / 4 .. 2 * pi]
    ]
    (eField protonOrigin)

simpleDipoleSodiumChloride :: ChargeDistribution
simpleDipoleSodiumChloride = simpleDipole (vec 0 0 2.992e-29) 2.36e-10

eFieldSodiumChloride :: VectorField
eFieldSodiumChloride = eField simpleDipoleSodiumChloride

eFieldPicSimpleDipole :: IO ()
eFieldPicSimpleDipole =
  vfGrad
    (** 0.2)
    (\(y, z) -> cart 0 (3e-10 * y) (3e-10 * z))
    (\v -> (yComp v, zComp v))
    "eFieldPicSimpleDipole.png"
    20
    eFieldSodiumChloride

eFieldIdealDipole ::
  Vec -> -- electric dipole moment
  VectorField -- electric field
eFieldIdealDipole p r =
  let k = 1 / (4 * pi * epsilon0) -- SI units
      rMag = magnitude (rVF r)
      rUnit = rVF r ^/ rMag
   in k *^ (1 / rMag ** 3) *^ (3 *^ (p <.> rUnit) *^ rUnit ^-^ p)

eFieldPicIdealDipole :: IO ()
eFieldPicIdealDipole =
  vfGrad
    (** 0.2)
    (\(y, z) -> cart 0 (3e-10 * y) (3e-10 * z))
    (\v -> (yComp v, zComp v))
    "eFieldPicIdealDipole.png"
    20
    (eFieldIdealDipole kHat)

lineDipoleSodiumChloride :: ChargeDistribution
lineDipoleSodiumChloride = lineDipole (vec 0 0 2.99e-29) 2.36e-10

eFieldLineDipole :: VectorField
eFieldLineDipole = eField lineDipoleSodiumChloride

eFieldDiskCap :: VectorField
eFieldDiskCap = eField $ diskCap 0.05 0.04 2e-8

eFieldPicDiskCap :: IO ()
eFieldPicDiskCap =
  vfGrad
    (** 0.2)
    (\(x, z) -> cart (0.1 * x) 0 (0.1 * z))
    (\v -> (xComp v, zComp v))
    "eFieldPicDiskCap.png"
    20
    eFieldDiskCap
