module Integration where

import Types  

type Integration = (R -> R) -- function
                  -> R      -- lower limit
                  -> R      -- upper limit
                  -> R      -- result

integral :: R ->  Integration
integral dt f a b
  = sum [f t * dt | t <- [a+dt/2 , a+3*dt/2 .. b - dt/2]]

type AntiDerivative = R         -- initial value
                    -> (R -> R) -- function
                    -> (R -> R) -- antiderivative of the function

antiDerivative :: R -> AntiDerivative
antiDerivative dt v0 a t = v0 + integral dt a 0 t            
