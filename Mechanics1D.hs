{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# OPTIONS_GHC -Wall #-}

module Mechanics1D where

import Graphics.Gnuplot.Simple
import Newton2 (fAir)
import Types (R)

type Time = R

type TimeStep = R

type Mass = R

type Position = R

type Velocity = R

type Force = R

type State1D = (Time, Position, Velocity)

newtonSecond1D ::
  Mass ->
  [State1D -> Force] -> -- force funcs
  State1D -> -- current state
  (R, R, R) -- deriv of state
newtonSecond1D m fs (t, x0, v0) =
  let fNet = sum [f (t, x0, v0) | f <- fs]
      acc = fNet / m
   in (1, v0, acc)

euler1D ::
  R -> -- time step dt
  (State1D -> (R, R, R)) -> -- differential equation
  (State1D -> State1D) --state update function
euler1D dt deriv (t0, x0, v0) =
  let (_, _, dvdt) = deriv (t0, x0, v0)
      t1 = t0 + dt
      x1 = x0 + v0 * dt
      v1 = v0 + dvdt * dt
   in (t1, x1, v1)

updateTXV ::
  R -> -- time interval dt
  Mass ->
  [State1D -> Force] -> -- list of force funcs
  (State1D -> State1D) -- state update func
updateTXV dt m fs = euler1D dt (newtonSecond1D m fs)

statesTXV ::
  R -> -- time interval dt
  Mass ->
  State1D -> -- initial state
  [State1D -> Force] -> -- list of force funcs
  [State1D] -- infinite list of states
statesTXV dt m txv0 fs = iterate (updateTXV dt m fs) txv0

-- assume that dt is the same between adjacent pairs
velocity1D ::
  [State1D] -> -- infinite list
  (Time -> Velocity) -- velocity function
velocity1D sts t =
  let (t0, _, _) = sts !! 0
      (t1, _, _) = sts !! 1
      dt = t1 - t0
      numSteps = abs $ round (t / dt)
      (_, _, v0) = sts !! numSteps
   in v0

velocityFtxv ::
  R -> -- time step
  Mass ->
  State1D -> -- initial state
  [State1D -> Force] -> -- list of force funcs
  (Time -> Velocity) -- velocity function
velocityFtxv dt m txv0 fs = velocity1D (statesTXV dt m txv0 fs)

-- assume that dt is the same between adjacent pairs
position1D ::
  [State1D] -> -- infinite list
  (Time -> Position) -- position function
position1D sts t =
  let (t0, _, _) = sts !! 0
      (t1, _, _) = sts !! 1
      dt = t1 - t0
      numSteps = abs $ round (t / dt)
      (_, x0, _) = sts !! numSteps
   in x0

positionFtxv ::
  R -> -- time step
  Mass ->
  State1D -> -- initial state
  [State1D -> Force] -> -- list of force funcs
  (Time -> Position) -- position function
positionFtxv dt m txv0 fs = position1D (statesTXV dt m txv0 fs)

springForce :: R -> State1D -> Force
springForce k (_, x0, _) = -k * x0

dampedHOForces :: [State1D -> Force]
dampedHOForces =
  [ springForce 0.8,
    \(_, _, v0) -> fAir 2 1.225 (pi * 0.02 ** 2) v0,
    \_ -> -0.0027 * 9.80665
  ]

dampedHOStates :: [State1D]
dampedHOStates = statesTXV 0.001 0.0027 (0.0, 0.1, 0.0) dampedHOForces

dampedHOGraph :: IO ()
dampedHOGraph =
  plotPath
    [ Title "Ping Pong Ball on a Slinky",
      XLabel "Time (s)",
      YLabel "Position (m)",
      PNG "dho.png",
      Key Nothing
    ]
    [(t, x) | (t, x, _) <- take 3000 dampedHOStates]

pingpongPosition :: Time -> Velocity
pingpongPosition = positionFtxv 0.001 0.0027 (0, 0.1, 0) dampedHOForces

dampedHOGraph2 :: IO ()
dampedHOGraph2 =
  plotFunc
    [ Title "Ping Pong Ball on a Slinky",
      XLabel "Time (s)",
      YLabel "Position (m)",
      PNG "dho1.png",
      Key Nothing
    ]
    [0, 0.01 .. 3]
    pingpongPosition

pingpongVelocity :: Time -> Velocity
pingpongVelocity = velocityFtxv 0.001 0.0027 (0, 0.1, 0) dampedHOForces

dampedHOGraph3 :: IO ()
dampedHOGraph3 =
  plotFunc
    [ Title "Ping Pong Ball on a Slinky",
      XLabel "Time (s)",
      YLabel "Velocity (m/s)",
      PNG "dho2.png",
      Key Nothing
    ]
    [0, 0.01 .. 3]
    pingpongVelocity

eulerCromer1D ::
  R -> -- time step dt
  (State1D -> (R, R, R)) -> -- differential equation
  (State1D -> State1D) -- state update function
eulerCromer1D dt deriv (t0, x0, v0) =
  let (_, _, dvdt) = deriv (t0, x0, v0)
      t1 = t0 + dt
      x1 = x0 + v1 * dt
      v1 = v0 + dvdt * dt
   in (t1, x1, v1)

updateTXVEC ::
  R -> -- time interval dt
  Mass ->
  [State1D -> Force] -> -- list of force funcs
  (State1D -> State1D) -- state update func
updateTXVEC dt m fs = eulerCromer1D dt (newtonSecond1D m fs)

type UpdateFunction s = s -> s

type DifferentialEquation s ds = s -> ds

type NumericalMethod s ds = DifferentialEquation s ds -> UpdateFunction s

solver :: NumericalMethod s ds -> DifferentialEquation s ds -> s -> [s]
solver method = iterate . method

class RealVectorSpace ds where
  (+++) :: ds -> ds -> ds
  scale :: R -> ds -> ds

instance RealVectorSpace (R, R, R) where
  (dtdt0, dxdt0, dvdt0) +++ (dtdt1, dxdt1, dvdt1) =
    (dtdt0 + dtdt1, dxdt0 + dxdt1, dvdt0 + dvdt1)
  scale w (dtdt0, dxdt0, dvdt0) = (w * dtdt0, w * dxdt0, w * dvdt0)

class RealVectorSpace ds => Diff s ds where
  shift :: R -> ds -> s -> s

instance Diff State1D (R, R, R) where
  shift dt (dtdt, dxdt, dvdt) (t, x, v) =
    (t + dtdt * dt, x + dxdt * dt, v + dvdt * dt)

euler :: Diff s ds => R -> (s -> ds) -> s -> s
euler dt deriv st0 = shift dt (deriv st0) st0

rungeKutta4 :: Diff s ds => R -> (s -> ds) -> s -> s
rungeKutta4 dt deriv st0 =
  let m0 = deriv st0
      m1 = deriv (shift (dt / 2) m0 st0)
      m2 = deriv (shift (dt / 2) m1 st0)
      m3 = deriv (shift dt m2 st0)
   in shift (dt / 6) (m0 +++ m1 +++ m1 +++ m2 +++ m2 +++ m3) st0

exponential :: DifferentialEquation (R, R, R) (R, R, R)
exponential (_, x0, v0) = (1, v0, x0)
