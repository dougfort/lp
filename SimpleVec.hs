{-# OPTIONS -Wall #-}

module SimpleVec where

type R = Double

infixl 6 ^+^

infixl 6 ^-^

infixr 7 *^

infixl 7 ^*

infixr 7 ^/

infixr 7 <.>

infixl 7 ><

data Vec = Vec
  { xComp :: R, -- x component
    yComp :: R, -- y component
    zComp :: R -- z xomponent
  }
  deriving (Eq)

instance Show Vec where
  show (Vec x y z) =
    "vec " ++ showDouble x ++ " "
      ++ showDouble y
      ++ " "
      ++ showDouble z

showDouble :: R -> String
showDouble x
  | x < 0 = "(" ++ show x ++ ")"
  | otherwise = show x

-- Form a vector by giving its x, y and z components:
vec ::
  R -> -- x component
  R -> -- y component
  R -> -- z component
  Vec
vec = Vec

iHat :: Vec
iHat = vec 1 0 0

jHat :: Vec
jHat = vec 0 1 0

kHat :: Vec
kHat = vec 0 0 1

zeroV :: Vec
zeroV = vec 0 0 0

negateV :: Vec -> Vec
negateV (Vec ax ay az) = Vec (-ax) (-ay) (-az)

(^+^) :: Vec -> Vec -> Vec
Vec ax ay az ^+^ Vec bx by bz = Vec (ax + bx) (ay + by) (az + bz)

(^-^) :: Vec -> Vec -> Vec
Vec ax ay az ^-^ Vec bx by bz = Vec (ax - bx) (ay - by) (az - bz)

sumV :: [Vec] -> Vec
sumV = foldr (^+^) zeroV

(*^) :: R -> Vec -> Vec
c *^ Vec ax ay az = Vec (c * ax) (c * ay) (c * az)

(^*) :: Vec -> R -> Vec
Vec ax ay az ^* c = Vec (c * ax) (c * ay) (c * az)

(<.>) :: Vec -> Vec -> R
Vec ax ay az <.> Vec bx by bz = ax * bx + ay * by + az * bz

(><) :: Vec -> Vec -> Vec
Vec ax ay az >< Vec bx by bz =
  Vec
    (ay * bz - az * by)
    (az * bx - ax * bz)
    (ax * by - ay * bx)

(^/) :: Vec -> R -> Vec
Vec ax ay az ^/ c = Vec (ax / c) (ay / c) (az / c)

magnitude :: Vec -> R
magnitude v = sqrt (v <.> v)

type PosVec = Vec
