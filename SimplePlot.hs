{-# OPTIONS -Wall #-}

module SimplePlot where


import Graphics.Gnuplot.Simple

type R = Double

tRange :: [R]
tRange = [0, 0.01..5]

yPos :: R -- y0
     -> R -- vy0
     -> R -- ay
     -> R -- t
     -> R -- y
yPos y0 vy0 ay t = y0 + vy0 * t + ay * t**2 / 2

plot1 :: IO ()
plot1 = plotFunc [Title "Projectile Motion"
                 ,XLabel "Time (S)"
                 ,YLabel "Height of projectile (m)"
                 ,PNG "projectie.png"
                 ,Key Nothing
                 ] tRange (yPos 0 20 (-9.8))
