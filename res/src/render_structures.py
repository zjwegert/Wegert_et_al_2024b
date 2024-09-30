import sys
# Paraview script to generate structures. This has been generated as a trace script
#  and modified to iterate over results and produce outputs.
# Author: ZJ Wegert
# Instructions: To execute from command line run `pvpython render_structures.py <RESULT PATH> <OUTPUT FOLDER NAME>`.
# Inputs:
RESULT_PATH = sys.argv[1] # E.g., '.../data/piezo'
OUTPUT_FOLDER = sys.argv[2] # E.g., 'data/structures/raw'

# Trace generated using paraview version 5.12.0-RC2
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 12

#### import the simple module from the paraview
from paraview.simple import *
from os import listdir, mkdir
from os.path import isfile, isdir, join, dirname, realpath
import re
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

OUTPUT_FOLDER = join(dirname(realpath(__file__)),'..',OUTPUT_FOLDER)
try:
  mkdir(OUTPUT_FOLDER)
except FileExistsError:
  pass
res_dirs = listdir(RESULT_PATH)
for res_dir in res_dirs:
  res_path=join(RESULT_PATH,res_dir)
  if not isdir(res_path):
    print("Skipping",res_dir,"as not a directory")
    continue
  dir_files = listdir(res_path)
  # Assume only one ptvu - replace with .vtu for serial structure
  ptvu_search = list(map(lambda x: re.search(r'out(\d+)\.([A-Za-z]+)',x),dir_files))
  ptvu_file = next(filter(lambda x: x!=None,ptvu_search),None)
  if ptvu_file == None:
    print("Error: could not find PTVU/TVU file in",res_dir)
    continue
  file_ext = ptvu_file[2]
  ptvu_file = ptvu_file.string

  #### PARAVIEW TRACE SCRIPT
  # create a new 'XML Partitioned Unstructured Grid Reader'
  # Replace with XMLUnstructuredGridReader for serial structure
  if file_ext == 'pvtu':
    outpvtu = XMLPartitionedUnstructuredGridReader(registrationName='outpvtu', FileName=[join(res_path,ptvu_file)])
  elif file_ext == 'vtu':
    outpvtu = XMLUnstructuredGridReader(registrationName='outpvtu', FileName=[join(res_path,ptvu_file)])
  else:
    print("Error: File extension",file_ext,"in",res_dir,"is not currently supported.")
    continue
  outpvtu.CellArrayStatus = ['piece', 'cell']
  outpvtu.PointArrayStatus = ['|∇(φ)|', 'φ', 'H(φ)']
  outpvtu.TimeArray = 'TimeValue'

  # Properties modified on outpvtu
  outpvtu.TimeArray = 'None'

  # get active view
  renderView1 = GetActiveViewOrCreate('RenderView')

  # show data in view
  outpvtuDisplay = Show(outpvtu, renderView1, 'UnstructuredGridRepresentation')

  # trace defaults for the display properties.
  outpvtuDisplay.Selection = None
  outpvtuDisplay.Representation = 'Surface'
  outpvtuDisplay.ColorArrayName = [None, '']
  outpvtuDisplay.LookupTable = None
  outpvtuDisplay.MapScalars = 1
  outpvtuDisplay.MultiComponentsMapping = 0
  outpvtuDisplay.InterpolateScalarsBeforeMapping = 1
  outpvtuDisplay.UseNanColorForMissingArrays = 0
  outpvtuDisplay.Opacity = 1.0
  outpvtuDisplay.PointSize = 2.0
  outpvtuDisplay.LineWidth = 1.0
  outpvtuDisplay.RenderLinesAsTubes = 0
  outpvtuDisplay.RenderPointsAsSpheres = 0
  outpvtuDisplay.Interpolation = 'Gouraud'
  outpvtuDisplay.Specular = 0.0
  outpvtuDisplay.SpecularColor = [1.0, 1.0, 1.0]
  outpvtuDisplay.SpecularPower = 100.0
  outpvtuDisplay.Luminosity = 0.0
  outpvtuDisplay.Ambient = 0.0
  outpvtuDisplay.Diffuse = 1.0
  outpvtuDisplay.Roughness = 0.3
  outpvtuDisplay.Metallic = 0.0
  outpvtuDisplay.EdgeTint = [1.0, 1.0, 1.0]
  outpvtuDisplay.Anisotropy = 0.0
  outpvtuDisplay.AnisotropyRotation = 0.0
  outpvtuDisplay.BaseIOR = 1.5
  outpvtuDisplay.CoatStrength = 0.0
  outpvtuDisplay.CoatIOR = 2.0
  outpvtuDisplay.CoatRoughness = 0.0
  outpvtuDisplay.CoatColor = [1.0, 1.0, 1.0]
  outpvtuDisplay.SelectTCoordArray = 'None'
  outpvtuDisplay.SelectNormalArray = 'None'
  outpvtuDisplay.SelectTangentArray = 'None'
  outpvtuDisplay.Texture = None
  outpvtuDisplay.RepeatTextures = 1
  outpvtuDisplay.InterpolateTextures = 0
  outpvtuDisplay.SeamlessU = 0
  outpvtuDisplay.SeamlessV = 0
  outpvtuDisplay.UseMipmapTextures = 0
  outpvtuDisplay.ShowTexturesOnBackface = 1
  outpvtuDisplay.BaseColorTexture = None
  outpvtuDisplay.NormalTexture = None
  outpvtuDisplay.NormalScale = 1.0
  outpvtuDisplay.CoatNormalTexture = None
  outpvtuDisplay.CoatNormalScale = 1.0
  outpvtuDisplay.MaterialTexture = None
  outpvtuDisplay.OcclusionStrength = 1.0
  outpvtuDisplay.AnisotropyTexture = None
  outpvtuDisplay.EmissiveTexture = None
  outpvtuDisplay.EmissiveFactor = [1.0, 1.0, 1.0]
  outpvtuDisplay.FlipTextures = 0
  outpvtuDisplay.EdgeOpacity = 1.0
  outpvtuDisplay.BackfaceRepresentation = 'Follow Frontface'
  outpvtuDisplay.BackfaceAmbientColor = [1.0, 1.0, 1.0]
  outpvtuDisplay.BackfaceOpacity = 1.0
  outpvtuDisplay.Position = [0.0, 0.0, 0.0]
  outpvtuDisplay.Scale = [1.0, 1.0, 1.0]
  outpvtuDisplay.Orientation = [0.0, 0.0, 0.0]
  outpvtuDisplay.Origin = [0.0, 0.0, 0.0]
  outpvtuDisplay.CoordinateShiftScaleMethod = 'Always Auto Shift Scale'
  outpvtuDisplay.Pickable = 1
  outpvtuDisplay.Triangulate = 0
  outpvtuDisplay.UseShaderReplacements = 0
  outpvtuDisplay.ShaderReplacements = ''
  outpvtuDisplay.NonlinearSubdivisionLevel = 1
  outpvtuDisplay.MatchBoundariesIgnoringCellOrder = 0
  outpvtuDisplay.UseDataPartitions = 0
  outpvtuDisplay.OSPRayUseScaleArray = 'All Approximate'
  outpvtuDisplay.OSPRayScaleArray = 'H(φ)'
  outpvtuDisplay.OSPRayScaleFunction = 'Piecewise Function'
  outpvtuDisplay.OSPRayMaterial = 'None'
  outpvtuDisplay.Assembly = ''
  outpvtuDisplay.BlockSelectors = ['/']
  outpvtuDisplay.BlockColors = []
  outpvtuDisplay.BlockOpacities = []
  outpvtuDisplay.Orient = 0
  outpvtuDisplay.OrientationMode = 'Direction'
  outpvtuDisplay.SelectOrientationVectors = 'None'
  outpvtuDisplay.Scaling = 0
  outpvtuDisplay.ScaleMode = 'No Data Scaling Off'
  outpvtuDisplay.ScaleFactor = 0.1
  outpvtuDisplay.SelectScaleArray = 'H(φ)'
  outpvtuDisplay.GlyphType = 'Arrow'
  outpvtuDisplay.UseGlyphTable = 0
  outpvtuDisplay.GlyphTableIndexArray = 'H(φ)'
  outpvtuDisplay.UseCompositeGlyphTable = 0
  outpvtuDisplay.UseGlyphCullingAndLOD = 0
  outpvtuDisplay.LODValues = []
  outpvtuDisplay.ColorByLODIndex = 0
  outpvtuDisplay.GaussianRadius = 0.005
  outpvtuDisplay.ShaderPreset = 'Sphere'
  outpvtuDisplay.CustomTriangleScale = 3
  outpvtuDisplay.CustomShader = """ // This custom shader code define a gaussian blur
  // Please take a look into vtkSMPointGaussianRepresentation.cxx
  // for other custom shader examples
  //VTK::Color::Impl
    float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
    float gaussian = exp(-0.5*dist2);
    opacity = opacity*gaussian;
  """
  outpvtuDisplay.Emissive = 0
  outpvtuDisplay.ScaleByArray = 0
  outpvtuDisplay.SetScaleArray = ['POINTS', 'H(φ)']
  outpvtuDisplay.ScaleArrayComponent = ''
  outpvtuDisplay.UseScaleFunction = 1
  outpvtuDisplay.ScaleTransferFunction = 'Piecewise Function'
  outpvtuDisplay.OpacityByArray = 0
  outpvtuDisplay.OpacityArray = ['POINTS', 'H(φ)']
  outpvtuDisplay.OpacityArrayComponent = ''
  outpvtuDisplay.OpacityTransferFunction = 'Piecewise Function'
  outpvtuDisplay.DataAxesGrid = 'Grid Axes Representation'
  outpvtuDisplay.SelectionCellLabelBold = 0
  outpvtuDisplay.SelectionCellLabelColor = [0.0, 1.0, 0.0]
  outpvtuDisplay.SelectionCellLabelFontFamily = 'Arial'
  outpvtuDisplay.SelectionCellLabelFontFile = ''
  outpvtuDisplay.SelectionCellLabelFontSize = 18
  outpvtuDisplay.SelectionCellLabelItalic = 0
  outpvtuDisplay.SelectionCellLabelJustification = 'Left'
  outpvtuDisplay.SelectionCellLabelOpacity = 1.0
  outpvtuDisplay.SelectionCellLabelShadow = 0
  outpvtuDisplay.SelectionPointLabelBold = 0
  outpvtuDisplay.SelectionPointLabelColor = [1.0, 1.0, 0.0]
  outpvtuDisplay.SelectionPointLabelFontFamily = 'Arial'
  outpvtuDisplay.SelectionPointLabelFontFile = ''
  outpvtuDisplay.SelectionPointLabelFontSize = 18
  outpvtuDisplay.SelectionPointLabelItalic = 0
  outpvtuDisplay.SelectionPointLabelJustification = 'Left'
  outpvtuDisplay.SelectionPointLabelOpacity = 1.0
  outpvtuDisplay.SelectionPointLabelShadow = 0
  outpvtuDisplay.PolarAxes = 'Polar Axes Representation'
  outpvtuDisplay.ScalarOpacityFunction = None
  outpvtuDisplay.ScalarOpacityUnitDistance = 0.017320508075688773
  outpvtuDisplay.UseSeparateOpacityArray = 0
  outpvtuDisplay.OpacityArrayName = ['POINTS', 'H(φ)']
  outpvtuDisplay.OpacityComponent = ''
  outpvtuDisplay.SelectMapper = 'Projected tetra'
  outpvtuDisplay.SamplingDimensions = [128, 128, 128]
  outpvtuDisplay.UseFloatingPointFrameBuffer = 1
  outpvtuDisplay.SelectInputVectors = [None, '']
  outpvtuDisplay.NumberOfSteps = 40
  outpvtuDisplay.StepSize = 0.25
  outpvtuDisplay.NormalizeVectors = 1
  outpvtuDisplay.EnhancedLIC = 1
  outpvtuDisplay.ColorMode = 'Blend'
  outpvtuDisplay.LICIntensity = 0.8
  outpvtuDisplay.MapModeBias = 0.0
  outpvtuDisplay.EnhanceContrast = 'Off'
  outpvtuDisplay.LowLICContrastEnhancementFactor = 0.0
  outpvtuDisplay.HighLICContrastEnhancementFactor = 0.0
  outpvtuDisplay.LowColorContrastEnhancementFactor = 0.0
  outpvtuDisplay.HighColorContrastEnhancementFactor = 0.0
  outpvtuDisplay.AntiAlias = 0
  outpvtuDisplay.MaskOnSurface = 1
  outpvtuDisplay.MaskThreshold = 0.0
  outpvtuDisplay.MaskIntensity = 0.0
  outpvtuDisplay.MaskColor = [0.5, 0.5, 0.5]
  outpvtuDisplay.GenerateNoiseTexture = 0
  outpvtuDisplay.NoiseType = 'Gaussian'
  outpvtuDisplay.NoiseTextureSize = 128
  outpvtuDisplay.NoiseGrainSize = 2
  outpvtuDisplay.MinNoiseValue = 0.0
  outpvtuDisplay.MaxNoiseValue = 0.8
  outpvtuDisplay.NumberOfNoiseLevels = 1024
  outpvtuDisplay.ImpulseNoiseProbability = 1.0
  outpvtuDisplay.ImpulseNoiseBackgroundValue = 0.0
  outpvtuDisplay.NoiseGeneratorSeed = 1
  outpvtuDisplay.CompositeStrategy = 'AUTO'
  outpvtuDisplay.UseLICForLOD = 0
  outpvtuDisplay.WriteLog = ''

  # init the 'Piecewise Function' selected for 'OSPRayScaleFunction'
  outpvtuDisplay.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
  outpvtuDisplay.OSPRayScaleFunction.UseLogScale = 0

  # init the 'Arrow' selected for 'GlyphType'
  outpvtuDisplay.GlyphType.TipResolution = 6
  outpvtuDisplay.GlyphType.TipRadius = 0.1
  outpvtuDisplay.GlyphType.TipLength = 0.35
  outpvtuDisplay.GlyphType.ShaftResolution = 6
  outpvtuDisplay.GlyphType.ShaftRadius = 0.03
  outpvtuDisplay.GlyphType.Invert = 0

  # init the 'Piecewise Function' selected for 'ScaleTransferFunction'
  outpvtuDisplay.ScaleTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
  outpvtuDisplay.ScaleTransferFunction.UseLogScale = 0

  # init the 'Piecewise Function' selected for 'OpacityTransferFunction'
  outpvtuDisplay.OpacityTransferFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
  outpvtuDisplay.OpacityTransferFunction.UseLogScale = 0

  # init the 'Grid Axes Representation' selected for 'DataAxesGrid'
  outpvtuDisplay.DataAxesGrid.XTitle = 'X Axis'
  outpvtuDisplay.DataAxesGrid.YTitle = 'Y Axis'
  outpvtuDisplay.DataAxesGrid.ZTitle = 'Z Axis'
  outpvtuDisplay.DataAxesGrid.XTitleFontFamily = 'Arial'
  outpvtuDisplay.DataAxesGrid.XTitleFontFile = ''
  outpvtuDisplay.DataAxesGrid.XTitleBold = 0
  outpvtuDisplay.DataAxesGrid.XTitleItalic = 0
  outpvtuDisplay.DataAxesGrid.XTitleFontSize = 12
  outpvtuDisplay.DataAxesGrid.XTitleShadow = 0
  outpvtuDisplay.DataAxesGrid.XTitleOpacity = 1.0
  outpvtuDisplay.DataAxesGrid.YTitleFontFamily = 'Arial'
  outpvtuDisplay.DataAxesGrid.YTitleFontFile = ''
  outpvtuDisplay.DataAxesGrid.YTitleBold = 0
  outpvtuDisplay.DataAxesGrid.YTitleItalic = 0
  outpvtuDisplay.DataAxesGrid.YTitleFontSize = 12
  outpvtuDisplay.DataAxesGrid.YTitleShadow = 0
  outpvtuDisplay.DataAxesGrid.YTitleOpacity = 1.0
  outpvtuDisplay.DataAxesGrid.ZTitleFontFamily = 'Arial'
  outpvtuDisplay.DataAxesGrid.ZTitleFontFile = ''
  outpvtuDisplay.DataAxesGrid.ZTitleBold = 0
  outpvtuDisplay.DataAxesGrid.ZTitleItalic = 0
  outpvtuDisplay.DataAxesGrid.ZTitleFontSize = 12
  outpvtuDisplay.DataAxesGrid.ZTitleShadow = 0
  outpvtuDisplay.DataAxesGrid.ZTitleOpacity = 1.0
  outpvtuDisplay.DataAxesGrid.FacesToRender = 63
  outpvtuDisplay.DataAxesGrid.CullBackface = 0
  outpvtuDisplay.DataAxesGrid.CullFrontface = 1
  outpvtuDisplay.DataAxesGrid.ShowGrid = 0
  outpvtuDisplay.DataAxesGrid.ShowEdges = 1
  outpvtuDisplay.DataAxesGrid.ShowTicks = 1
  outpvtuDisplay.DataAxesGrid.LabelUniqueEdgesOnly = 1
  outpvtuDisplay.DataAxesGrid.AxesToLabel = 63
  outpvtuDisplay.DataAxesGrid.XLabelFontFamily = 'Arial'
  outpvtuDisplay.DataAxesGrid.XLabelFontFile = ''
  outpvtuDisplay.DataAxesGrid.XLabelBold = 0
  outpvtuDisplay.DataAxesGrid.XLabelItalic = 0
  outpvtuDisplay.DataAxesGrid.XLabelFontSize = 12
  outpvtuDisplay.DataAxesGrid.XLabelShadow = 0
  outpvtuDisplay.DataAxesGrid.XLabelOpacity = 1.0
  outpvtuDisplay.DataAxesGrid.YLabelFontFamily = 'Arial'
  outpvtuDisplay.DataAxesGrid.YLabelFontFile = ''
  outpvtuDisplay.DataAxesGrid.YLabelBold = 0
  outpvtuDisplay.DataAxesGrid.YLabelItalic = 0
  outpvtuDisplay.DataAxesGrid.YLabelFontSize = 12
  outpvtuDisplay.DataAxesGrid.YLabelShadow = 0
  outpvtuDisplay.DataAxesGrid.YLabelOpacity = 1.0
  outpvtuDisplay.DataAxesGrid.ZLabelFontFamily = 'Arial'
  outpvtuDisplay.DataAxesGrid.ZLabelFontFile = ''
  outpvtuDisplay.DataAxesGrid.ZLabelBold = 0
  outpvtuDisplay.DataAxesGrid.ZLabelItalic = 0
  outpvtuDisplay.DataAxesGrid.ZLabelFontSize = 12
  outpvtuDisplay.DataAxesGrid.ZLabelShadow = 0
  outpvtuDisplay.DataAxesGrid.ZLabelOpacity = 1.0
  outpvtuDisplay.DataAxesGrid.XAxisNotation = 'Mixed'
  outpvtuDisplay.DataAxesGrid.XAxisPrecision = 2
  outpvtuDisplay.DataAxesGrid.XAxisUseCustomLabels = 0
  outpvtuDisplay.DataAxesGrid.XAxisLabels = []
  outpvtuDisplay.DataAxesGrid.YAxisNotation = 'Mixed'
  outpvtuDisplay.DataAxesGrid.YAxisPrecision = 2
  outpvtuDisplay.DataAxesGrid.YAxisUseCustomLabels = 0
  outpvtuDisplay.DataAxesGrid.YAxisLabels = []
  outpvtuDisplay.DataAxesGrid.ZAxisNotation = 'Mixed'
  outpvtuDisplay.DataAxesGrid.ZAxisPrecision = 2
  outpvtuDisplay.DataAxesGrid.ZAxisUseCustomLabels = 0
  outpvtuDisplay.DataAxesGrid.ZAxisLabels = []
  outpvtuDisplay.DataAxesGrid.UseCustomBounds = 0
  outpvtuDisplay.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

  # init the 'Polar Axes Representation' selected for 'PolarAxes'
  outpvtuDisplay.PolarAxes.Visibility = 0
  outpvtuDisplay.PolarAxes.Translation = [0.0, 0.0, 0.0]
  outpvtuDisplay.PolarAxes.Scale = [1.0, 1.0, 1.0]
  outpvtuDisplay.PolarAxes.Orientation = [0.0, 0.0, 0.0]
  outpvtuDisplay.PolarAxes.EnableCustomBounds = [0, 0, 0]
  outpvtuDisplay.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
  outpvtuDisplay.PolarAxes.EnableCustomRange = 0
  outpvtuDisplay.PolarAxes.CustomRange = [0.0, 1.0]
  outpvtuDisplay.PolarAxes.AutoPole = 1
  outpvtuDisplay.PolarAxes.PolarAxisVisibility = 1
  outpvtuDisplay.PolarAxes.RadialAxesVisibility = 1
  outpvtuDisplay.PolarAxes.DrawRadialGridlines = 1
  outpvtuDisplay.PolarAxes.PolarArcsVisibility = 1
  outpvtuDisplay.PolarAxes.DrawPolarArcsGridlines = 1
  outpvtuDisplay.PolarAxes.NumberOfRadialAxes = 0
  outpvtuDisplay.PolarAxes.DeltaAngleRadialAxes = 45.0
  outpvtuDisplay.PolarAxes.NumberOfPolarAxes = 5
  outpvtuDisplay.PolarAxes.DeltaRangePolarAxes = 0.0
  outpvtuDisplay.PolarAxes.CustomMinRadius = 1
  outpvtuDisplay.PolarAxes.MinimumRadius = 0.0
  outpvtuDisplay.PolarAxes.CustomAngles = 1
  outpvtuDisplay.PolarAxes.MinimumAngle = 0.0
  outpvtuDisplay.PolarAxes.MaximumAngle = 90.0
  outpvtuDisplay.PolarAxes.RadialAxesOriginToPolarAxis = 1
  outpvtuDisplay.PolarAxes.PolarArcResolutionPerDegree = 0.2
  outpvtuDisplay.PolarAxes.Ratio = 1.0
  outpvtuDisplay.PolarAxes.EnableOverallColor = 1
  outpvtuDisplay.PolarAxes.OverallColor = [1.0, 1.0, 1.0]
  outpvtuDisplay.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
  outpvtuDisplay.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
  outpvtuDisplay.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
  outpvtuDisplay.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
  outpvtuDisplay.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
  outpvtuDisplay.PolarAxes.PolarAxisTitleVisibility = 1
  outpvtuDisplay.PolarAxes.PolarAxisTitle = 'Radial Distance'
  outpvtuDisplay.PolarAxes.PolarAxisTitleLocation = 'Bottom'
  outpvtuDisplay.PolarAxes.PolarTitleOffset = [20.0, 20.0]
  outpvtuDisplay.PolarAxes.PolarLabelVisibility = 1
  outpvtuDisplay.PolarAxes.PolarLabelFormat = '%-#6.3g'
  outpvtuDisplay.PolarAxes.PolarLabelExponentLocation = 'Labels'
  outpvtuDisplay.PolarAxes.PolarLabelOffset = 10.0
  outpvtuDisplay.PolarAxes.PolarExponentOffset = 5.0
  outpvtuDisplay.PolarAxes.RadialTitleVisibility = 1
  outpvtuDisplay.PolarAxes.RadialTitleFormat = '%-#3.1f'
  outpvtuDisplay.PolarAxes.RadialTitleLocation = 'Bottom'
  outpvtuDisplay.PolarAxes.RadialTitleOffset = [20.0, 0.0]
  outpvtuDisplay.PolarAxes.RadialUnitsVisibility = 1
  outpvtuDisplay.PolarAxes.ScreenSize = 10.0
  outpvtuDisplay.PolarAxes.PolarAxisTitleOpacity = 1.0
  outpvtuDisplay.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
  outpvtuDisplay.PolarAxes.PolarAxisTitleFontFile = ''
  outpvtuDisplay.PolarAxes.PolarAxisTitleBold = 0
  outpvtuDisplay.PolarAxes.PolarAxisTitleItalic = 0
  outpvtuDisplay.PolarAxes.PolarAxisTitleShadow = 0
  outpvtuDisplay.PolarAxes.PolarAxisTitleFontSize = 12
  outpvtuDisplay.PolarAxes.PolarAxisLabelOpacity = 1.0
  outpvtuDisplay.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
  outpvtuDisplay.PolarAxes.PolarAxisLabelFontFile = ''
  outpvtuDisplay.PolarAxes.PolarAxisLabelBold = 0
  outpvtuDisplay.PolarAxes.PolarAxisLabelItalic = 0
  outpvtuDisplay.PolarAxes.PolarAxisLabelShadow = 0
  outpvtuDisplay.PolarAxes.PolarAxisLabelFontSize = 12
  outpvtuDisplay.PolarAxes.LastRadialAxisTextOpacity = 1.0
  outpvtuDisplay.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
  outpvtuDisplay.PolarAxes.LastRadialAxisTextFontFile = ''
  outpvtuDisplay.PolarAxes.LastRadialAxisTextBold = 0
  outpvtuDisplay.PolarAxes.LastRadialAxisTextItalic = 0
  outpvtuDisplay.PolarAxes.LastRadialAxisTextShadow = 0
  outpvtuDisplay.PolarAxes.LastRadialAxisTextFontSize = 12
  outpvtuDisplay.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
  outpvtuDisplay.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
  outpvtuDisplay.PolarAxes.SecondaryRadialAxesTextFontFile = ''
  outpvtuDisplay.PolarAxes.SecondaryRadialAxesTextBold = 0
  outpvtuDisplay.PolarAxes.SecondaryRadialAxesTextItalic = 0
  outpvtuDisplay.PolarAxes.SecondaryRadialAxesTextShadow = 0
  outpvtuDisplay.PolarAxes.SecondaryRadialAxesTextFontSize = 12
  outpvtuDisplay.PolarAxes.EnableDistanceLOD = 1
  outpvtuDisplay.PolarAxes.DistanceLODThreshold = 0.7
  outpvtuDisplay.PolarAxes.EnableViewAngleLOD = 1
  outpvtuDisplay.PolarAxes.ViewAngleLODThreshold = 0.7
  outpvtuDisplay.PolarAxes.SmallestVisiblePolarAngle = 0.5
  outpvtuDisplay.PolarAxes.PolarTicksVisibility = 1
  outpvtuDisplay.PolarAxes.ArcTicksOriginToPolarAxis = 1
  outpvtuDisplay.PolarAxes.TickLocation = 'Both'
  outpvtuDisplay.PolarAxes.AxisTickVisibility = 1
  outpvtuDisplay.PolarAxes.AxisMinorTickVisibility = 0
  outpvtuDisplay.PolarAxes.AxisTickMatchesPolarAxes = 1
  outpvtuDisplay.PolarAxes.DeltaRangeMajor = 1.0
  outpvtuDisplay.PolarAxes.DeltaRangeMinor = 0.5
  outpvtuDisplay.PolarAxes.ArcTickVisibility = 1
  outpvtuDisplay.PolarAxes.ArcMinorTickVisibility = 0
  outpvtuDisplay.PolarAxes.ArcTickMatchesRadialAxes = 1
  outpvtuDisplay.PolarAxes.DeltaAngleMajor = 10.0
  outpvtuDisplay.PolarAxes.DeltaAngleMinor = 5.0
  outpvtuDisplay.PolarAxes.TickRatioRadiusSize = 0.02
  outpvtuDisplay.PolarAxes.PolarAxisMajorTickSize = 0.0
  outpvtuDisplay.PolarAxes.PolarAxisTickRatioSize = 0.3
  outpvtuDisplay.PolarAxes.PolarAxisMajorTickThickness = 1.0
  outpvtuDisplay.PolarAxes.PolarAxisTickRatioThickness = 0.5
  outpvtuDisplay.PolarAxes.LastRadialAxisMajorTickSize = 0.0
  outpvtuDisplay.PolarAxes.LastRadialAxisTickRatioSize = 0.3
  outpvtuDisplay.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
  outpvtuDisplay.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
  outpvtuDisplay.PolarAxes.ArcMajorTickSize = 0.0
  outpvtuDisplay.PolarAxes.ArcTickRatioSize = 0.3
  outpvtuDisplay.PolarAxes.ArcMajorTickThickness = 1.0
  outpvtuDisplay.PolarAxes.ArcTickRatioThickness = 0.5
  outpvtuDisplay.PolarAxes.Use2DMode = 0
  outpvtuDisplay.PolarAxes.UseLogAxis = 0

  # reset view to fit data
  renderView1.ResetCamera(False, 0.9)

  # get the material library
  materialLibrary1 = GetMaterialLibrary()

  # update the view to ensure updated data information
  renderView1.Update()
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
  renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
  renderView1.CameraParallelScale = 0.8660254037844386
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
  renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
  renderView1.CameraParallelScale = 0.8660254037844386
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
  renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
  renderView1.CameraParallelScale = 0.8660254037844386
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
  renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
  renderView1.CameraParallelScale = 0.8660254037844386
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
  renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
  renderView1.CameraParallelScale = 0.8660254037844386
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
  renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
  renderView1.CameraParallelScale = 0.8660254037844386
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
  renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
  renderView1.CameraParallelScale = 0.8660254037844386
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
  renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
  renderView1.CameraParallelScale = 0.8660254037844386
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [0.5, 0.5, 3.8460652149512318]
  renderView1.CameraFocalPoint = [0.5, 0.5, 0.5]
  renderView1.CameraParallelScale = 0.8660254037844386

  # Properties modified on renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Hide orientation axes
  renderView1.OrientationAxesVisibility = 0
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # create a new 'Iso Volume'
  isoVolume1 = IsoVolume(registrationName='IsoVolume1', Input=outpvtu)
  isoVolume1.InputScalars = ['POINTS', 'φ']
  isoVolume1.ThresholdRange = [-10.0, 0.00]

  # show data in view
  isoVolume1Display = Show(isoVolume1, renderView1, 'UnstructuredGridRepresentation')

  # hide data in view
  Hide(outpvtu, renderView1)

  # show color bar/color legend
  isoVolume1Display.SetScalarBarVisibility(renderView1, True)

  # update the view to ensure updated data information
  renderView1.Update()
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # hide color bar/color legend
  isoVolume1Display.SetScalarBarVisibility(renderView1, False)
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # turn off scalar coloring
  ColorBy(isoVolume1Display, None)
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.EnableRayTracing = 1
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.LightScale = 0.75
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.AmbientSamples = 10
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.SamplesPerPixel = 10
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.ProgressivePasses = 10
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # get layout
  layout1 = GetLayout()

  # layout/tab size in pixels
  layout1.SetSize(806, 794)

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # save screenshot
  SaveScreenshot(join(OUTPUT_FOLDER,res_dir+'_STRUCTURE.png'), renderView1, 16, ImageResolution=[1612, 1588],
      FontScaling='Scale fonts proportionally',
      OverrideColorPalette='',
      StereoMode='No change',
      TransparentBackground=0,
      SaveInBackground=0,
      EmbedParaViewState=0,
      # PNG options
      CompressionLevel='5',
      MetaData=['Application', 'ParaView'])
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  #================================================================
  # addendum: following script captures some of the application
  # state to faithfully reproduce the visualization during playback
  #================================================================

  #--------------------------------
  # saving layout sizes for layouts

  # layout/tab size in pixels
  layout1.SetSize(806, 794)

  #-----------------------------------
  # saving camera placements for views

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Render void phase
  isoVolume1.ThresholdRange = [0.00, 10.0]

  # show data in view
  isoVolume1Display = Show(isoVolume1, renderView1, 'UnstructuredGridRepresentation')

  # hide data in view
  Hide(outpvtu, renderView1)

  # show color bar/color legend
  isoVolume1Display.SetScalarBarVisibility(renderView1, True)

  # update the view to ensure updated data information
  renderView1.Update()
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # hide color bar/color legend
  isoVolume1Display.SetScalarBarVisibility(renderView1, False)
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # turn off scalar coloring
  ColorBy(isoVolume1Display, None)
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.EnableRayTracing = 1
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.LightScale = 0.75
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.AmbientSamples = 10
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.SamplesPerPixel = 10
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # Properties modified on renderView1
  renderView1.ProgressivePasses = 10
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # get layout
  layout1 = GetLayout()

  # layout/tab size in pixels
  layout1.SetSize(806, 794)

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  # save screenshot
  SaveScreenshot(join(OUTPUT_FOLDER,res_dir+'_STRUCTURE_VOID.png'), renderView1, 16, ImageResolution=[1612, 1588],
      FontScaling='Scale fonts proportionally',
      OverrideColorPalette='',
      StereoMode='No change',
      TransparentBackground=0,
      SaveInBackground=0,
      EmbedParaViewState=0,
      # PNG options
      CompressionLevel='5',
      MetaData=['Application', 'ParaView'])
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168
  # Adjust camera

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  #================================================================
  # addendum: following script captures some of the application
  # state to faithfully reproduce the visualization during playback
  #================================================================

  #--------------------------------
  # saving layout sizes for layouts

  # layout/tab size in pixels
  layout1.SetSize(806, 794)

  #-----------------------------------
  # saving camera placements for views

  # current camera placement for renderView1
  renderView1.CameraPosition = [1.8012884763430275, -1.855070382959957, 1.3548173321272308]
  renderView1.CameraFocalPoint = [0.4939566564121575, 0.450950497718811, 0.39085793323172985]
  renderView1.CameraViewUp = [-0.17943003074334082, 0.29109636940833566, 0.9397168551136718]
  renderView1.CameraParallelScale = 0.7872958216222168

  Disconnect()
  Connect()

# ##--------------------------------------------
# ## You may need to add some code at the end of this python script depending on your usage, eg:
# #
# ## Render all views to see them appears
# # RenderAllViews()
# #
# ## Interact with the view, usefull when running from pvpython
# # Interact()
# #
# ## Save a screenshot of the active view
# # SaveScreenshot("path/to/screenshot.png")
# #
# ## Save a screenshot of a layout (multiple splitted view)
# # SaveScreenshot("path/to/screenshot.png", GetLayout())
# #
# ## Save all "Extractors" from the pipeline browser
# # SaveExtracts()
# #
# ## Save a animation of the current active view
# # SaveAnimation()
# #
# ## Please refer to the documentation of paraview.simple
# ## https://kitware.github.io/paraview-docs/latest/python/paraview.simple.html
# ##--------------------------------------------