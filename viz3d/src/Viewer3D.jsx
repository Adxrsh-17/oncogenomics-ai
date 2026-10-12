import { useEffect, useRef, useState, useCallback } from 'react';

const STYLES = { cartoon: 'cartoon', surface: 'surface', 'ball-stick': 'ball-stick' };

export default function Viewer3D({ pdbStr, accentColor, mutResidue, mode, style, showDrug, showMutation, spinning }) {
  const containerRef = useRef(null);
  const viewerRef    = useRef(null);
  const spinRef      = useRef(spinning);
  const [loading, setLoading]   = useState(true);
  const [error,   setError]     = useState(null);

  const buildScene = useCallback(() => {
    const viewer = viewerRef.current;
    if (!viewer || !pdbStr) return;

    viewer.removeAllModels();
    viewer.removeAllSurfaces();
    viewer.removeAllShapes();
    viewer.clear();

    try {
      viewer.addModel(pdbStr, 'pdb');

      if (style === 'cartoon') {
        viewer.setStyle({}, { cartoon: { color: 'spectrum', thickness: 0.4 } });
      } else if (style === 'surface') {
        viewer.setStyle({}, { cartoon: { color: 'spectrum', opacity: 0.25, thickness: 0.3 } });
        viewer.addSurface(window.$3Dmol.SurfaceType.VDW, { opacity: 0.65, colorscheme: 'whiteCarbon' });
      } else {
        viewer.setStyle({}, { stick: { radius: 0.14, colorscheme: 'Jmol' }, sphere: { scale: 0.22 } });
      }

      if (showMutation && mutResidue) {
        viewer.setStyle(
          { resi: mutResidue },
          {
            sphere: { radius: 1.6, color: accentColor, opacity: 0.9 },
            stick:  { radius: 0.3,  color: accentColor },
          }
        );
      }

      if (mode === 'crystal' && showDrug) {
        viewer.setStyle(
          { hetflag: true },
          { stick: { radius: 0.28, colorscheme: 'Jmol' }, sphere: { scale: 0.42, colorscheme: 'Jmol' } }
        );
      }

      viewer.zoomTo();
      viewer.render();
      setLoading(false);

      if (spinRef.current) viewer.spin('y', 0.4);
      else viewer.spin(false);
    } catch (e) {
      setError(e.message);
      setLoading(false);
    }
  }, [pdbStr, accentColor, mutResidue, mode, style, showDrug, showMutation]);

  // Init viewer once
  useEffect(() => {
    if (!containerRef.current || viewerRef.current) return;
    const wait3dmol = setInterval(() => {
      if (!window.$3Dmol) return;
      clearInterval(wait3dmol);
      const el = containerRef.current;
      if (!el) return;
      viewerRef.current = window.$3Dmol.createViewer(el, {
        backgroundColor: '#030508',
        antialias: true,
        cartoonQuality: 6,
      });
      buildScene();
    }, 50);
    return () => clearInterval(wait3dmol);
  }, []); // eslint-disable-line

  // Rebuild on scene params change
  useEffect(() => {
    if (viewerRef.current) {
      setLoading(true);
      buildScene();
    }
  }, [buildScene]);

  // Spin toggle without full rebuild
  useEffect(() => {
    spinRef.current = spinning;
    if (!viewerRef.current) return;
    if (spinning) viewerRef.current.spin('y', 0.4);
    else viewerRef.current.spin(false);
  }, [spinning]);

  return (
    <div style={{ position: 'relative', width: '100%', height: '100%', background: '#030508', borderRadius: 12, overflow: 'hidden' }}>
      <div ref={containerRef} style={{ width: '100%', height: '100%' }} />
      {loading && !error && (
        <div style={overlayStyle}>
          <div style={spinnerStyle} />
          <div style={{ fontFamily: 'monospace', fontSize: 12, color: '#00d4ff', marginTop: 8 }}>
            Rendering structure…
          </div>
        </div>
      )}
      {error && (
        <div style={{ ...overlayStyle, color: '#ff6b8a', fontSize: 13, padding: 24, textAlign: 'center' }}>
          3Dmol error: {error}
        </div>
      )}
    </div>
  );
}

const overlayStyle = {
  position: 'absolute', inset: 0, display: 'flex', flexDirection: 'column',
  alignItems: 'center', justifyContent: 'center', background: '#030508',
};
const spinnerStyle = {
  width: 36, height: 36, border: '3px solid #1e3050', borderTopColor: '#00d4ff',
  borderRadius: '50%', animation: 'spin 0.8s linear infinite',
};
