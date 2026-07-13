const STEPS = [
  { label: 'Setup',      id: 'setup' },
  { label: 'Metacells',  id: 'metacells' },
  { label: 'SetDatExpr', id: 'datexpr' },
  { label: 'Soft Power', id: 'softpower', breakPoint: true },
  { label: 'Network',    id: 'network' },
  { label: 'Modules',    id: 'modules' },
]

export default function PipelineStepper({ activeStep }) {
  return (
    <nav className="stepper" aria-label="Pipeline stages">
      {STEPS.map((step, i) => {
        const state =
          i < activeStep  ? 'done'   :
          i === activeStep ? 'active' :
          activeStep < 3 && i > activeStep ? 'locked' : ''

        const cls = ['stepper-step', state, step.breakPoint ? 'break-point' : '']
          .filter(Boolean).join(' ')

        return (
          <span key={step.id} style={{ display: 'contents' }}>
            <span className={cls} aria-current={i === activeStep ? 'step' : undefined}>
              <span className="step-circle">
                {i < activeStep ? '✓' : i + 1}
              </span>
              {step.label}
            </span>
            {i < STEPS.length - 1 && (
              <span className="stepper-arrow" aria-hidden="true">›</span>
            )}
          </span>
        )
      })}
    </nav>
  )
}
