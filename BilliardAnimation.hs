{-# OPTIONS -Wall #-}

import MOExamples (billiardInitial, billiardPicture, billiardUpdate)
import Mechanics3D (simulateGloss)
import MultipleObjects (eulerCromerMPS)

main :: IO ()
main =
  simulateGloss
    1
    100
    billiardInitial
    billiardPicture
    (billiardUpdate eulerCromerMPS 30)
