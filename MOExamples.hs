{-# OPTIONS -Wall #-}

module MOExamples where

-- orient

import qualified Graphics.Gloss as G
import Graphics.Gnuplot.Simple
import Mechanics1D
  ( NumericalMethod,
    TimeStep,
    euler,
    rungeKutta4,
  )
import Mechanics3D
  ( HasTime (..),
    ParticleState (..),
    customLabel,
    defaultParticleState,
    disk,
    earthSurfaceGravity,
  )
import MultipleObjects
  ( DMultiParticleState (..),
    Force (..),
    MultiParticleState (..),
    TwoBodyForce,
    billiardForce,
    eulerCromerMPS,
    fixedLinearSpring,
    linearSpring,
    newtonSecondMPS,
    statesMPS,
    updateMPS,
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
    (*^),
    (^+^),
    (^-^),
  )

twoSpringForces :: [Force]
twoSpringForces =
  [ ExternalForce 0 (fixedLinearSpring 100 0.5 zeroV),
    InternalForce 0 1 (linearSpring 100 0.5),
    ExternalForce 0 earthSurfaceGravity,
    ExternalForce 1 earthSurfaceGravity
  ]

twoSpringInitial :: MultiParticleState
twoSpringInitial =
  MPS
    [ defaultParticleState
        { mass = 2,
          posVec = 0.4 *^ jHat ^-^ 0.3 *^ kHat
        },
      defaultParticleState
        { mass = 3,
          posVec = 0.4 *^ jHat ^-^ 0.8 *^ kHat
        }
    ]

twoSpringUpdate ::
  TimeStep ->
  MultiParticleState -> -- old state
  MultiParticleState -- new state
twoSpringUpdate dt = updateMPS (eulerCromerMPS dt) twoSpringForces

kineticEnergy :: ParticleState -> R
kineticEnergy st =
  let m = mass st
      v = magnitude (velocity st)
   in (1 / 2) * m * v ** 2

systemKE :: MultiParticleState -> R
systemKE (MPS sts) = sum [kineticEnergy st | st <- sts]

linearSpringPE ::
  R -> -- spring constant
  R -> -- equilibrium length
  ParticleState -> -- state of particle at one end of string
  ParticleState -> -- state of particle at other end of string
  R -- potential energy of the spring
linearSpringPE k re st1 st2 =
  let r1 = posVec st1
      r2 = posVec st2
      r21 = r2 ^-^ r1
      r21mag = magnitude r21
   in k * (r21mag - re) ** 2 / 2

-- z direction is towrd the sky
-- assumes SI units
earthSurfaceGravityPE :: ParticleState -> R
earthSurfaceGravityPE st =
  let g = 9.80665 -- m/s^2
      m = mass st
      z = zComp (posVec st)
   in m * g * z

twoSpringsPE :: MultiParticleState -> R
twoSpringsPE (MPS sts) =
  linearSpringPE 100 0.5 defaultParticleState (head sts)
    + linearSpringPE 100 0.5 (head sts) (sts !! 1)
    + earthSurfaceGravityPE (head sts)
    + earthSurfaceGravityPE (sts !! 1)

twoSpringsME :: MultiParticleState -> R
twoSpringsME mpst = systemKE mpst + twoSpringsPE mpst

ballRadius :: R
ballRadius = 0.03 -- 6 cm diameter = 0.03 m radius

billiardForces :: R -> [Force]
billiardForces k = [InternalForce 0 1 (billiardForce k (2 * ballRadius))]

billiardDiffEq :: R -> MultiParticleState -> DMultiParticleState
billiardDiffEq k = newtonSecondMPS $ billiardForces k

billiardUpdate ::
  (TimeStep -> NumericalMethod MultiParticleState DMultiParticleState) ->
  R -> -- k
  TimeStep -> -- dt
  (MultiParticleState -> MultiParticleState)
billiardUpdate nMethod k dt = updateMPS (nMethod dt) (billiardForces k)

billiardEvolver ::
  (TimeStep -> NumericalMethod MultiParticleState DMultiParticleState) ->
  R -> -- k
  TimeStep -> -- dt
  (MultiParticleState -> [MultiParticleState])
billiardEvolver nMethod k dt = statesMPS (nMethod dt) (billiardForces k)

billiardInitial :: MultiParticleState
billiardInitial =
  let ballMass = 0.160 -- 160 g
   in MPS
        [ defaultParticleState
            { mass = ballMass,
              posVec = zeroV,
              velocity = 0.2 *^ iHat
            },
          defaultParticleState
            { mass = ballMass,
              posVec = iHat ^+^ 0.02 *^ jHat,
              velocity = zeroV
            }
        ]

billiardStates ::
  (TimeStep -> NumericalMethod MultiParticleState DMultiParticleState) ->
  R -> -- k
  TimeStep -> -- dt
  [MultiParticleState]
billiardStates nMethod k dt =
  statesMPS (nMethod dt) (billiardForces k) billiardInitial

billiardStatesFinite ::
  (TimeStep -> NumericalMethod MultiParticleState DMultiParticleState) ->
  R -> -- k
  TimeStep -> -- dt
  [MultiParticleState]
billiardStatesFinite nMethod k dt =
  takeWhile (\st -> timeOf st <= 10) (billiardStates nMethod k dt)

momentum :: ParticleState -> Vec
momentum st =
  let m = mass st
      v = velocity st
   in m *^ v

systemP :: MultiParticleState -> Vec
systemP (MPS sts) = sumV [momentum st | st <- sts]

percentChangePMag :: [MultiParticleState] -> R
percentChangePMag mpsts =
  let p0 = systemP (head mpsts)
      p1 = systemP (last mpsts)
   in 100 * magnitude (p1 ^-^ p0) / magnitude p0

sigFigs :: Int -> R -> Float
sigFigs n x =
  let expon :: Int
      expon = floor (logBase 10 x) - n + 1
      toInt :: R -> Int
      toInt = round
   in (10 ^^ expon *) $ fromIntegral $ toInt (10 ^^ (-expon) * x)

data Justification = LJ | RJ deriving (Show)

data Table a = Table Justification [[a]]

instance Show a => Show (Table a) where
  show (Table j xss) =
    let pairWithLength x = let str = show x in (str, length str)
        pairss = map (map pairWithLength) xss
        maxLength = maximum (map (maximum . map snd) pairss)
        showPair (str, len) =
          case j of
            LJ -> str ++ replicate (maxLength + 1 - len) ' '
            RJ -> replicate (maxLength + 1 - len) ' ' ++ str
        showLine pairs = concatMap showPair pairs ++ "\n"
     in init $ concatMap showLine pairss

pTable ::
  (TimeStep -> NumericalMethod MultiParticleState DMultiParticleState) ->
  [R] -> -- ks
  [TimeStep] -> -- dts
  Table Float
pTable nMethod ks dts =
  Table
    LJ
    [ [ sigFigs 2 $ percentChangePMag (billiardStatesFinite nMethod k dt)
        | dt <- dts
      ]
      | k <- ks
    ]

pTableEu ::
  [R] -> -- ks
  [TimeStep] -> -- dts
  Table Float
pTableEu = pTable euler
