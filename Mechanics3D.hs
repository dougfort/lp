{-# OPTIONS -Wall #-}
{-# LANGUAGE MultiParamTypeClasses #-}

module Mechanics3D where

import qualified Graphics.Gloss as G
import Graphics.Gnuplot.Simple
  ( Aspect (..),
    Attribute (..),
    plotFunc,
    plotPaths,
  )
import Mechanics1D
  ( Diff (..),
    NumericalMethod,
    RealVectorSpace (..),
    Time,
    TimeStep,
    rungeKutta4,
    solver,
  )
import SimpleVec
  ( PosVec,
    R,
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
    (<.>),
    (><),
    (^*),
    (^+^),
    (^-^),
    (^/),
  )
import SpatialMath
  ( Euler (..),
    V3 (..)
  )
import qualified Vis as V

data ParticleState = ParticleState
  { mass :: R,
    charge :: R,
    time :: R,
    posVec :: Vec,
    velocity :: Vec
  }
  deriving (Show)

defaultParticleState :: ParticleState
defaultParticleState =
  ParticleState
    { mass = 1,
      charge = 1,
      time = 0,
      posVec = zeroV,
      velocity = zeroV
    }

type OneBodyForce = ParticleState -> Vec

data DParticleState = DParticleState
  { dmdt :: R,
    dqdt :: R,
    dtdt :: R,
    drdt :: Vec,
    dvdt :: Vec
  }
  deriving (Show)

instance RealVectorSpace DParticleState where
  dst1 +++ dst2 =
    DParticleState
      { dmdt = dmdt dst1 + dmdt dst2,
        dqdt = dqdt dst1 + dqdt dst2,
        dtdt = dtdt dst1 + dtdt dst2,
        drdt = drdt dst1 ^+^ drdt dst2,
        dvdt = dvdt dst1 ^+^ dvdt dst2
      }
  scale w dst =
    DParticleState
      { dmdt = w * dmdt dst,
        dqdt = w * dqdt dst,
        dtdt = w * dtdt dst,
        drdt = w *^ drdt dst,
        dvdt = w *^ dvdt dst
      }

instance Diff ParticleState DParticleState where
  shift dt dps (ParticleState m q t r v) =
    ParticleState
      (m + dmdt dps * dt)
      (q + dqdt dps * dt)
      (t + dtdt dps * dt)
      (r ^+^ drdt dps ^* dt)
      (v ^+^ dvdt dps ^* dt)

class HasTime s where
  timeOf :: s -> Time

instance HasTime ParticleState where
  timeOf = time

newtonSecondPS :: [OneBodyForce] -> ParticleState -> DParticleState -- difeq
newtonSecondPS fs st =
  let fNet = sumV [f st | f <- fs]
      m = mass st
      v = velocity st
      acc = fNet ^/ m
   in DParticleState
        { dmdt = 0,
          dqdt = 0,
          dtdt = 1,
          drdt = v,
          dvdt = acc
        }

-- z direction is toward the sky
-- assumes SI units
earthSurfaceGravity :: OneBodyForce
earthSurfaceGravity st =
  let g = 9.80665 -- m/s^2
   in (-mass st * g) *^ kHat

-- origin at center of sun
-- assumes SI units
sunGravity :: OneBodyForce
sunGravity (ParticleState m _q _t r _v) =
  let bigG = 6.67408e-11 -- N m^2/kg^2
      sunMass = 1.98848e30 -- kg
   in (-bigG * sunMass * m) *^ r ^/ magnitude r ** 3

airResistance ::
  R -> -- drag coefficient
  R -> -- air density
  R -> -- cross sectional area of object
  OneBodyForce
airResistance drag rho area (ParticleState _m _q _t _r v) =
  (-0.5 * drag * rho * area * magnitude v) *^ v

windForce ::
  Vec -> -- wind velocity
  R -> -- drag coefficient
  R -> -- air density
  R -> -- cross sectional area of object
  OneBodyForce
windForce vWind drag rho area (ParticleState _m _q _t _r v) =
  let vRel = v ^-^ vWind
   in (-0.5 * drag * rho * area * magnitude vRel) *^ vRel

uniformLorentzForce ::
  Vec -> -- E
  Vec -> -- B
  OneBodyForce
uniformLorentzForce vE vB (ParticleState _m q _t _r v) =
  q *^ (vE ^+^ v >< vB)

eulerCromerPS ::
  TimeStep -> -- dt for stepping
  NumericalMethod ParticleState DParticleState
eulerCromerPS dt deriv st =
  let t = time st
      r = posVec st
      v = velocity st
      dst = deriv st
      acc = dvdt dst
      v' = v ^+^ acc ^* dt
   in st
        { time = t + dt,
          posVec = r ^+^ v' ^* dt,
          velocity = v ^+^ acc ^* dt
        }

statesPS ::
  NumericalMethod ParticleState DParticleState ->
  [OneBodyForce] ->
  (ParticleState -> [ParticleState])
statesPS method = iterate . method . newtonSecondPS

updatePS ::
  NumericalMethod ParticleState DParticleState ->
  [OneBodyForce] ->
  (ParticleState -> ParticleState)
updatePS method = method . newtonSecondPS

positionPS ::
  NumericalMethod ParticleState DParticleState ->
  [OneBodyForce] -> -- list of force funcs
  ParticleState -> -- initial state
  (Time -> PosVec) -- position function
positionPS method fs st t =
  let states = statesPS method fs st
      dt = time (states !! 1) - time (head states)
      numSteps = abs $ round (t / dt)
      st1 = solver method (newtonSecondPS fs) st !! numSteps
   in posVec st1

simulateGloss ::
  R -> -- time scale factor
  Int -> -- animation rate
  s -> -- initial state
  (s -> G.Picture) ->
  (TimeStep -> s -> s) ->
  IO ()
simulateGloss tsFactor rate initialState picFunc updateFunc =
  G.simulate
    (G.InWindow "" (1000, 750) (10, 10))
    G.black
    rate
    initialState
    picFunc
    (\_ -> updateFunc . (* tsFactor) . realToFrac)

simulateVis ::
  HasTime s =>
  R -> -- time scale factor
  Int -> -- animation rate
  s -> -- initial state
  (s -> V.VisObject R) ->
  (TimeStep -> s -> s) ->
  IO ()
simulateVis tsFactor rate initialState picFunc updateFunc =
  let visUpdateFunc ta st = let dtp = tsFactor * realToFrac ta - timeOf st in updateFunc dtp st
   in V.simulate
        V.defaultOpts
        (1 / fromIntegral rate)
        initialState
        (orient . picFunc)
        visUpdateFunc

v3FromVec :: Vec -> V3 R
v3FromVec v = V3 x y z 
  where 
    x = xComp v
    y = yComp v
    z = zComp v

orient :: V.VisObject R -> V.VisObject R
orient = V.RotEulerDeg (Euler 270 180 0)

halleyUpdate ::
  TimeStep ->
  (ParticleState -> ParticleState)
halleyUpdate dt =
  updatePS (eulerCromerPS dt) [sunGravity]

halleyInitial :: ParticleState
halleyInitial =
  ParticleState
    { mass = 2.2e14,
      charge = 0,
      time = 0,
      posVec = 8.766e10 *^ iHat,
      velocity = 54569 *^ jHat
    }

disk :: Float -> G.Picture
disk radius = G.ThickCircle (radius / 2) radius

baseballForces :: [OneBodyForce]
baseballForces =
  let area = pi * (0.074 / 2) ** 2
   in [earthSurfaceGravity, airResistance 0.3 1.225 area]

baseballTrajectory ::
  R -> --time step
  R -> -- initial speed
  R -> -- launch angle in degrees
  [(R, R)] -- (y, z) pairs
baseballTrajectory dt v0 thetaDeg =
  let thetaRad = thetaDeg * pi / 180
      vy0 = v0 * cos thetaRad
      vz0 = v0 * sin thetaRad
      initialState =
        ParticleState
          { mass = 0.145,
            charge = 0,
            time = 0,
            posVec = zeroV,
            velocity = vec 0 vy0 vz0
          }
   in trajectory $
        zGEO $
          statesPS (eulerCromerPS dt) baseballForces initialState

zGEO :: [ParticleState] -> [ParticleState]
zGEO = takeWhile (\(ParticleState _ _ _ r _) -> zComp r >= 0)

trajectory :: [ParticleState] -> [(R, R)]
trajectory sts = [(yComp r, zComp r) | ParticleState _ _ _ r _ <- sts]

baseballRange ::
  R -> -- time step
  R -> -- initial speed
  R -> -- launch angle in degrees
  R -- range
baseballRange dt v0 thetaDeg =
  let (y, _) = last $ baseballTrajectory dt v0 thetaDeg
   in y

baseballRangeGraph :: IO ()
baseballRangeGraph =
  plotFunc
    [ Title "Range for baseball hit at 45 m/s",
      XLabel "Angle above horizontal (degrees)",
      YLabel "Horizontal range (m)",
      PNG "baseballrange.png",
      Key Nothing
    ]
    [10, 11 .. 80]
    $ baseballRange 0.01 45

bestAngle :: (R, R)
bestAngle =
  maximum
    [ (baseballRange 0.01 45 thetaDeg, thetaDeg)
      | thetaDeg <- [30, 31 .. 60]
    ]

projectileUpdate ::
  TimeStep ->
  ParticleState -> -- old state
  ParticleState -- new state
projectileUpdate dt =
  updatePS (eulerCromerPS dt) baseballForces

projectileInitial :: [String] -> ParticleState
projectileInitial [] = error "Please supply initial speed and angle"
projectileInitial [_] = error "Please supply initial speed and angle"
projectileInitial (_ : _ : _ : _) = error "First argument is speed. Second argument is degrees."
projectileInitial (arg1 : arg2 : _) =
  let v0 = read arg1 :: R -- initial speed m/s
      angleDeg = read arg2 :: R -- initial angle degrees
      theta = angleDeg * pi / 180 -- in radians
   in defaultParticleState
        { mass = 0.145,
          charge = 0,
          time = 0,
          posVec = zeroV,
          velocity = vec 0 (v0 * cos theta) (v0 * sin theta)
        }

relativityPS ::
  [OneBodyForce] ->
  (ParticleState -> DParticleState) -- a differential equation
relativityPS fs st =
  let fNet = sumV [f st | f <- fs]
      c = 299792458 -- m/s
      m = mass st
      v = velocity st
      u = v ^/ c
      acc = sqrt (1 - u <.> u) *^ (fNet ^-^ (fNet <.> u) *^ u) ^/ m
   in DParticleState
        { dmdt = 0,
          dqdt = 0,
          dtdt = 1,
          drdt = v,
          dvdt = acc
        }

constantForcePlot :: IO ()
constantForcePlot =
  let year = 365.25 * 24 * 60 * 60 -- seconds
      c = 299792458 -- m/s
      method = rungeKutta4 1000
      forces = [const (10 *^ iHat)]
      initialState = defaultParticleState {mass = 1}
      newtonStates = solver method (newtonSecondPS forces) initialState
      relativityStates = solver method (relativityPS forces) initialState
      newtonTVs =
        [ (time st / year, xComp (velocity st) / c)
          | st <- takeWhile tle1yr newtonStates
        ]
      relativityTVs =
        [ (time st / year, xComp (velocity st) / c)
          | st <- takeWhile tle1yr relativityStates
        ]
   in plotPaths
        [ Key Nothing,
          Title "Response to a constant force",
          XLabel "Time (years)",
          YLabel "Velocity (multiples of x)",
          PNG "constantforce.PNG",
          customLabel (0.1, 1) "mass = 1 kg",
          customLabel (0.1, 0.9) "force = 10 N",
          customLabel (0.5, 0.7) "Newtonian",
          customLabel (0.8, 0.6) "relativistic"
        ]
        [newtonTVs, relativityTVs]

customLabel :: (R, R) -> String -> Attribute
customLabel (x, y) label =
  Custom "label" ["\"" ++ label ++ "\"" ++ " at " ++ show x ++ "," ++ show y]

tle1yr :: ParticleState -> Bool
tle1yr st =
  let year = 365.25 * 24 * 60 * 60
   in time st <= year

circularPlot :: IO ()
circularPlot =
  let c = 299792458 -- m/s
      method = rungeKutta4 1e-9
      forces = [uniformLorentzForce zeroV kHat] -- 1 T
      initialState =
        defaultParticleState
          { mass = 1.672621898e-27,
            charge = 1.602176634e-19,
            velocity = 0.8 *^ c *^ jHat
          }
      newtonStates = solver method (newtonSecondPS forces) initialState
      relativityStates = solver method (relativityPS forces) initialState
      newtonXYs = [(xComp (posVec st), yComp (posVec st)) | st <- take 100 newtonStates]
      relativityXYs = [(xComp (posVec st), yComp (posVec st)) | st <- take 120 relativityStates]
   in plotPaths
        [ Key Nothing,
          Aspect (Ratio 1),
          Title "Proton in a 1-T magnetic field",
          XLabel "x (m)",
          YLabel "y (m)",
          PNG "circularComp.png",
          customLabel (0.5, 4.5) "v = 0.8 c",
          customLabel (2.5, 0.0) "Newtonian",
          customLabel (3.0, 3.5) "relativistic"
        ]
        [newtonXYs, relativityXYs]
