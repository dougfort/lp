{-# OPTIONS -Wall #-}

import Diagrams.Prelude (black, yellow)
import Maxwell (initialStateFDTD, jGaussian, makeEpng, stateUpdate)

main :: IO ()
main =
  let dt = 0.02e-9 -- 0.02 ns time step
      numTimeSteps = 719
   in mapM_
        (makeEpng (yellow, black))
        ( zip [0 .. numTimeSteps] $
            iterate (stateUpdate dt jGaussian) (initialStateFDTD 0.108)
        )
