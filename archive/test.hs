{-# LANGUAGE FlexibleInstances #-}

import Data.Map
import Control.Monad.State



data WQConfig   = WQConfig { wqname :: String
                           , wqport :: Int
                           } deriving (Show)

data WorkQueue  = WorkQueue                    deriving Show

data Topology   = Topology                     deriving Show


data Coord      = XYZ (Double, Double, Double) deriving Show

data Structure  = Structure [Coord]            deriving Show

data Molecule   = Molecule Topology Structure  deriving Show

type CellID     = Int
type Weight     = Double
type Color      = String

data Cells      = Cells { centers :: [Structure]
                        , colors  :: [Color]
                        } deriving Show

data Walker     = Walker { wStructure :: Structure 
                         , wWeight    :: Weight
                         , wColor     :: Color
                         , wCell      :: CellID
                         } deriving Show

data Transition = TransitionMatrix deriving Show

data Resample   = Resample { rWalkers     :: [Walker]
                           , rTransitions :: Transition
                           } deriving Show

-- keep the type checker happy
instance Show (Resample -> IO Resample) where show _ = "<resample function>"


data AWE        = AWE { aweWqconfig :: WQConfig
                      , aweCells    :: Cells
                      , aweWalkers  :: [Walker]
                      , aweTopology :: Topology
                      , aweRasample :: Resample -> IO Resample
                      } deriving Show



resample :: Resample -> IO Resample
resample r = do
  -- update weights and colors
  return r